module neuralSurface_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill
  use trainedMLP_class,   only : trainedMLP

  implicit none
  private

  character(*), parameter :: TYPE_NAME = 'neuralSurface'

  !!
  !! Neural SDF surface
  !!
  !! An implicit surface defined by a trained MLP that approximates a signed
  !! distance function F(r):
  !!   F(r) < 0  =>  r is inside the surface  (negative halfspace)
  !!   F(r) > 0  =>  r is outside the surface (positive halfspace)
  !!
  !! Intended for delta (Woodcock) tracking only: distance() returns INF, so it
  !! is not valid for surface-tracking modes. going() uses a single forward
  !! finite-difference step along the particle direction to resolve the
  !! halfspace on the surface boundary.
  !!
  !! Memory: the MLP weight arrays are allocatable components of the trainedMLP
  !! and are freed by kill() and, as a safety net, by the finaliser. init()
  !! kills any existing MLP before loading, so re-initialisation is safe.
  !!
  !! See misclassClerk_class (Tallies/TallyClerks) for a halfspace
  !! misclassification diagnostic against a reference region.
  !!
  !! Private Members:
  !!   mlp       -> Trained MLP loaded from the weight file at init time
  !!   geomScale -> Geometric scale factor: physical coordinates are divided by
  !!                this before MLP evaluation, letting a unit-sphere weight
  !!                file represent a sphere of arbitrary radius. Default 1.0.
  !!
  !! Interface:
  !!   myType      -> Return surface type name
  !!   init        -> Build from a dictionary (reads the weight file)
  !!   boundingBox -> MLP training box scaled by geomScale
  !!   evaluate    -> Signed surface value at a point
  !!   distance    -> Stub (INF); delta-tracking-only surface
  !!   going       -> Finite-difference halfspace resolution on the boundary
  !!   kill        -> Free the MLP and reset to the uninitialised state
  !!
  !! Sample Dictionary Input:
  !!   ns { type neuralSurface;
  !!        id 1;
  !!        weightFile "sphere_weights.bin";
  !!        # geometricScale 1.0; #
  !!      }
  !!
  type, public, extends(surface) :: neuralSurface
    private
    type(trainedMLP) :: mlp
    real(defReal)    :: geomScale = ONE
  contains
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
    final     :: finaliseNeuralSurface
  end type neuralSurface

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for details
  !!
  pure function myType(self) result(str)
    class(neuralSurface), intent(in) :: self
    character(:), allocatable        :: str

    str = TYPE_NAME

  end function myType

  !!
  !! Initialise from a dictionary
  !!
  !! Reads 'id' and 'weightFile', loads the MLP from that binary or text weight
  !! file, and reads the optional 'geometricScale'. Any MLP already loaded (e.g.
  !! after kill()) is released first.
  !!
  !! See surface_inter for details
  !!
  !! Errors:
  !!   fatalError if id < 1, geometricScale <= 0, or the weight file cannot be read
  !!
  subroutine init(self, dict)
    class(neuralSurface), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: id
    character(pathLen)                  :: weightFile
    character(100), parameter :: Here = 'init (neuralSurface_class.f90)'

    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Invalid surface id: '//numToChar(id))

    call dict % get(weightFile, 'weightFile')

    ! Release any existing MLP weight arrays before (re-)loading
    if (self % mlp % isInit) call self % mlp % kill()

    call self % mlp % load(trim(weightFile))
    call self % setId(id)

    ! Optional geometric scale: divides physical coordinates before MLP
    ! evaluation, letting a unit-sphere weight file represent any radius.
    call dict % getOrDefault(self % geomScale, 'geometricScale', ONE)
    if (self % geomScale <= ZERO) call fatalError(Here, 'geometricScale must be positive')

  end subroutine init

  !!
  !! Return the axis-aligned bounding box
  !!
  !! The MLP training bounding box scaled by geomScale.
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(neuralSurface), intent(in) :: self
    real(defReal), dimension(6)      :: aabb

    aabb(1:3) = self % mlp % bboxMin * self % geomScale
    aabb(4:6) = self % mlp % bboxMax * self % geomScale

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! Delegates to the pure MLP forward pass, scaling world coordinates by
  !! 1 / geomScale on the way in and the result by geomScale on the way out.
  !! Negative return = inside (negative halfspace).
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(neuralSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    c = self % mlp % evaluate(r / self % geomScale) * self % geomScale

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! Stub: always INF. neuralSurface is for delta (Woodcock) tracking only;
  !! an analytic distance is not available for an arbitrary MLP.
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(neuralSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d

    d = INF

  end function distance

  !!
  !! Return TRUE if the particle is going into the +ve halfspace
  !!
  !! Uses a single forward finite-difference step along u: the particle is
  !! moving toward the +ve halfspace if F(r + FD_STEP * u) > 0.
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(hs)
    class(neuralSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    real(defReal), parameter :: FD_STEP = 1.0e-7_defReal

    hs = self % evaluate(r + FD_STEP * u) > ZERO

  end function going

  !!
  !! Return to the uninitialised state and free the MLP weight arrays
  !!
  !! See surface_inter for details
  !!
  elemental subroutine kill(self)
    class(neuralSurface), intent(inout) :: self

    call kill_super(self)
    call self % mlp % kill()
    self % geomScale = ONE

  end subroutine kill

  !!
  !! Finaliser: a safety net that frees the MLP if kill() was never called
  !!
  !! Called automatically by Fortran when a neuralSurface is deallocated or
  !! goes out of scope.
  !!
  subroutine finaliseNeuralSurface(self)
    type(neuralSurface), intent(inout) :: self

    call self % kill()

  end subroutine finaliseNeuralSurface

end module neuralSurface_class
