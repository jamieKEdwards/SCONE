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

  ! Module-level diagnostic counters and file handle
  integer(longInt), public              :: neuralSurf_nCalls    = 0_longInt
  integer(longInt), public              :: neuralSurf_nMisclass = 0_longInt
  integer(shortInt), parameter, private :: DIAG_UNIT = 97
  logical(defBool),  save,      private :: diagFileOpen = .false.

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
  !!   mlp        -> Trained MLP loaded from weight file at init time
  !!   geomScale  -> Geometric scale factor: physical coords are divided by this
  !!                 before MLP evaluation, allowing a unit-sphere weight file to
  !!                 represent a sphere of arbitrary radius. Default 1.0 (no scaling).
  !!   diagRefR2  -> Reference sphere R^2 for misclassification diagnostics;
  !!                 0 means diagnostic is disabled (default)
  !!
  type, public, extends(surface) :: neuralSurface
    private
    type(trainedMLP) :: mlp
    real(defReal)    :: geomScale  = ONE
    real(defReal)    :: diagRefR2  = ZERO
  contains
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: halfspace
    final     :: finaliseNeuralSurface
  end type neuralSurface

  public :: printNeuralDiagnostics

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
    real(defReal)                 :: diagRadius
    character(100), parameter :: Here = 'init (neuralSurface_class.f90)'

    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Invalid surface id: '//numToChar(id))

    call dict % get(weightFile, 'weightFile')

    ! Release any existing MLP weight arrays before (re-)loading
    if (self % mlp % isInit) call self % mlp % kill()

    call readMLPWeights(self % mlp, trim(weightFile))
    call self % setId(id)

    ! Optional geometric scale: divides physical coords before MLP evaluation.
    ! Allows a unit-sphere weight file to represent a sphere of arbitrary radius.
    if (dict % isPresent('geometricScale')) then
      call dict % get(self % geomScale, 'geometricScale')
      if (self % geomScale <= ZERO) &
        call fatalError(Here, 'geometricScale must be positive')
    else
      self % geomScale = ONE
    end if

    ! Optional diagnostic: reference sphere radius for misclassification counting
    if (dict % isPresent('diagRadius')) then
      call dict % get(diagRadius, 'diagRadius')
      self % diagRefR2 = diagRadius * diagRadius
      if (.not. diagFileOpen) then
        open(unit=DIAG_UNIT, file='neural_misclass.dat', status='replace', action='write')
        write(DIAG_UNIT, '(A)') '# x  y  z  neural_outside  sphere_outside  sdf'
        diagFileOpen = .true.
      end if
    else
      self % diagRefR2 = ZERO
    end if

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

    aabb(1:3) = self % mlp % bboxMin * self % geomScale
    aabb(4:6) = self % mlp % bboxMax * self % geomScale

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

    c = self % mlp % evaluate(r / self % geomScale) * self % geomScale

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

    hs = self % mlp % evaluate((r + FD_STEP * u) / self % geomScale) > ZERO

  end function going

  !!
  !! Return true if particle is in +ve halfspace
  !!
  !! Overrides the default surface_inter implementation to add runtime
  !! misclassification diagnostics when diagRadius is set in the input.
  !! For each call, the neural halfspace is compared against the analytic
  !! sphere and the module-level OMP ATOMIC counters are updated.
  !!
  !! See surface_inter for details
  !!
  function halfspace(self, r, u) result(hs)
    class(neuralSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    logical(defBool)                        :: sphere_hs
    real(defReal)                           :: c

    ! Evaluate neural SDF and determine halfspace
    c = self % mlp % evaluate(r / self % geomScale) * self % geomScale
    if (abs(c) < self % surfTol()) then
      hs = self % going(r, u)
    else
      hs = c > ZERO
    end if

    ! Misclassification diagnostic: compare with analytic sphere (if enabled)
    if (self % diagRefR2 > ZERO) then
      sphere_hs = (r(1)*r(1) + r(2)*r(2) + r(3)*r(3) - self % diagRefR2) > ZERO
      !$omp atomic
      neuralSurf_nCalls = neuralSurf_nCalls + 1_longInt
      if (hs .neqv. sphere_hs) then
        !$omp atomic
        neuralSurf_nMisclass = neuralSurf_nMisclass + 1_longInt
        if (diagFileOpen) then
          !$omp critical(neuralDiag)
          write(DIAG_UNIT, '(3ES16.8, 2L3, ES16.8)') r(1), r(2), r(3), hs, sphere_hs, c
          !$omp end critical(neuralDiag)
        end if
      end if
    end if

  end function halfspace

  !!
  !! Print misclassification diagnostic summary
  !!
  !! Prints total halfspace call count and misclassification count/rate.
  !! Does nothing if no calls were tallied (diagnostic not enabled or not used).
  !!
  subroutine printNeuralDiagnostics()
    real(defReal) :: pct

    if (neuralSurf_nCalls == 0_longInt) return

    pct = 100.0_defReal * real(neuralSurf_nMisclass, defReal) / real(neuralSurf_nCalls, defReal)

    print '(A)', ''
    print '(A)', '--- Neural Surface Halfspace Diagnostic ---'
    print '(A,I0)', '  Total halfspace calls : ', neuralSurf_nCalls
    print '(A,I0)', '  Misclassified calls   : ', neuralSurf_nMisclass
    print '(A,F8.4,A)', '  Misclassification rate: ', pct, ' %'
    if (diagFileOpen) then
      print '(A)', '  Details written to: neural_misclass.dat'
      write(DIAG_UNIT, '(A)')       '# ---- Summary ----'
      write(DIAG_UNIT, '(A,I0)')    '# Total calls:   ', neuralSurf_nCalls
      write(DIAG_UNIT, '(A,I0)')    '# Misclassified: ', neuralSurf_nMisclass
      write(DIAG_UNIT, '(A,F8.4,A)') '# Rate:          ', pct, ' %'
      close(DIAG_UNIT)
      diagFileOpen = .false.
    end if
    print '(A)', '-------------------------------------------'

  end subroutine printNeuralDiagnostics

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
