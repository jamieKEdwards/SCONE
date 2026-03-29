module neuralSurface_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface
  use mlpInference_mod,   only : trainedMLP
  use mlpWeightIO_mod,    only : readMLPWeights

  implicit none
  private

  character(*), parameter :: TYPE_NAME = 'neuralSurface'

  !!
  !! Neural SDF surface
  !!
  !! Represents an implicit surface defined by a trained MLP (Multi-Layer Perceptron).
  !! The MLP approximates a signed distance function F(r):
  !!   F(r) < 0  =>  r is inside the surface  (negative halfspace)
  !!   F(r) > 0  =>  r is outside the surface (positive halfspace)
  !!
  !! Intended for use with delta (Woodcock) tracking only.
  !! distance() returns INF — not valid for surface-tracking modes.
  !!
  !! going() uses a single forward finite-difference step along the
  !! particle direction to resolve the halfspace on the surface boundary.
  !!
  !! Memory management: weights are freed by the FINAL subroutine when the
  !! object is deallocated. init() will kill any existing MLP before loading
  !! new weights, so re-initialisation is safe after kill().
  !!
  !! Sample dictionary input:
  !!   ns { type neuralSurface;
  !!        id 1;
  !!        weightFile "sphere_weights.bin";
  !!      }
  !!
  !! Private Members:
  !!   mlp -> Trained MLP loaded from weight file at init time
  !!
  type, public, extends(surface) :: neuralSurface
    private
    type(trainedMLP) :: mlp
  contains
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
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
  !! Initialise from dictionary
  !!
  !! Reads 'id' and 'weightFile' entries from the dictionary.
  !! Loads the MLP from the binary or text weight file at 'weightFile'.
  !! If an MLP is already loaded (e.g. after kill()), it is released first.
  !!
  !! See surface_inter for details
  !!
  !! Errors:
  !!   fatalError if id < 1 or weight file cannot be read
  !!
  subroutine init(self, dict)
    class(neuralSurface), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)             :: id
    character(pathLen)            :: weightFile
    character(100), parameter :: Here = 'init (neuralSurface_class.f90)'

    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Invalid surface id: '//numToChar(id))

    call dict % get(weightFile, 'weightFile')

    ! Release any existing MLP weight arrays before (re-)loading
    if (self % mlp % isInit) call self % mlp % kill()

    call readMLPWeights(self % mlp, trim(weightFile))
    call self % setId(id)

  end subroutine init

  !!
  !! Return axis-aligned bounding box
  !!
  !! Returns the MLP training bounding box as the surface bbox.
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(neuralSurface), intent(in) :: self
    real(defReal), dimension(6)      :: aabb

    aabb(1:3) = self % mlp % bboxMin
    aabb(4:6) = self % mlp % bboxMax

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! Delegates to the pure MLP forward pass.
  !! Negative return = inside (negative halfspace).
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(neuralSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    c = self % mlp % evaluate(r)

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! Stub: always returns INF.
  !! neuralSurface is for delta (Woodcock) tracking only;
  !! analytic distance is not available for an arbitrary MLP.
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
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! Uses a single forward finite-difference step along u.
  !! If F(r + FD_STEP * u) > 0 the particle is moving toward the +ve halfspace.
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(hs)
    class(neuralSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    real(defReal), parameter :: FD_STEP = 1.0e-7_defReal

    hs = self % mlp % evaluate(r + FD_STEP * u) > ZERO

  end function going

  !!
  !! Finaliser: release MLP weight arrays when object is destroyed
  !!
  !! Called automatically by Fortran when a neuralSurface object is
  !! deallocated or goes out of scope. Ensures weight arrays are freed
  !! even if kill() was not explicitly called.
  !!
  subroutine finaliseNeuralSurface(self)
    type(neuralSurface), intent(inout) :: self

    if (self % mlp % isInit) call self % mlp % kill()

  end subroutine finaliseNeuralSurface

end module neuralSurface_class
