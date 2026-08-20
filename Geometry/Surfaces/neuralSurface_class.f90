module neuralSurface_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,        only : surface
  use mlpInference_mod,     only : trainedMLP
  use mlpWeightIO_mod,      only : readMLPWeights
  use sphere_class,         only : sphere        !! !!TO BE REMOVED!!
  use bezierShape_class,    only : bezierShape   !! !!TO BE REMOVED!!
  use bezierTwist_class,    only : bezierTwist   !! !!TO BE REMOVED!!
  use bezierVolume_class,   only : bezierVolume  !! !!TO BE REMOVED!! (diagnostic refSurfaces union only)

  implicit none
  private

  character(*), parameter :: TYPE_NAME = 'neuralSurface'

  !! !!TO BE REMOVED!! Module-level diagnostic counters and file handle
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
  !! distance() returns INF -- not valid for surface-tracking modes.
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
  !!
  type, public, extends(surface) :: neuralSurface
    private
    type(trainedMLP)           :: mlp
    real(defReal)              :: geomScale   = ONE
    !! !!TO BE REMOVED!! Misclassification diagnostic members. refSurfBody/
    !! Handle/Spout: union-of-3 reference for complex multi-shape geometry
    !! (the teapot), same convention as deepLSSurface_class.f90's
    !! identical addition -- see that module for the full derivation.
    real(defReal)              :: diagRefR2   = ZERO
    class(surface), pointer    :: refSurf       => null()
    class(surface), pointer    :: refSurfBody   => null()
    class(surface), pointer    :: refSurfHandle => null()
    class(surface), pointer    :: refSurfSpout  => null()
    logical(defBool)           :: refFlip     = .false.
    logical(defBool)           :: diagEnabled = .true.
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

  public :: printNeuralDiagnostics !! !!TO BE REMOVED!!

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
    !! !!TO BE REMOVED!! Diagnostic local variables
    real(defReal)                 :: diagRadius
    character(pathLen)            :: diagFile
    integer(shortInt)             :: refFlipInt
    integer(shortInt)             :: diagEnabledInt
    type(dictionary)              :: refSurfDict
    character(nameLen)            :: refType
    type(dictionary)              :: refSurfacesDict
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

    !! !!TO BE REMOVED!! diagFile default
    if (dict % isPresent('diagFile')) then
      call dict % get(diagFile, 'diagFile')
    else
      diagFile = 'neural_misclass.dat'
    end if

    !! !!TO BE REMOVED!! Optional diagnostic: reference sphere radius for misclassification counting
    if (dict % isPresent('diagRadius')) then
      call dict % get(diagRadius, 'diagRadius')
      self % diagRefR2 = diagRadius * diagRadius
      if (.not. diagFileOpen) then
        open(unit=DIAG_UNIT, file=trim(diagFile), status='replace', action='write')
        write(DIAG_UNIT, '(A)') '# x  y  z  neural_outside  sphere_outside  sdf'
        diagFileOpen = .true.
      end if
    else
      self % diagRefR2 = ZERO
    end if

    !! !!TO BE REMOVED!! Optional diagnostic: reference surface for misclassification comparison.
    !! Supports sphere, bezierShape, bezierTwist as refSurface types.
    if (dict % isPresent('refSurface')) then
      if (dict % isPresent('diagEnabled')) then
        call dict % get(diagEnabledInt, 'diagEnabled')
        self % diagEnabled = (diagEnabledInt /= 0)
      else
        self % diagEnabled = .true.
      end if
      if (associated(self % refSurf)) then
        call self % refSurf % kill()
        deallocate(self % refSurf)
      end if
      call dict % get(refSurfDict, 'refSurface')
      call refSurfDict % get(refType, 'type')
      select case (trim(refType))
        case ('sphere')
          allocate(sphere :: self % refSurf)
        case ('bezierShape')
          allocate(bezierShape :: self % refSurf)
        case ('bezierTwist')
          allocate(bezierTwist :: self % refSurf)
        case default
          call fatalError(Here, 'Unsupported refSurface type: '//trim(refType))
      end select
      call self % refSurf % init(refSurfDict)
      if (dict % isPresent('refFlip')) then
        call dict % get(refFlipInt, 'refFlip')
        self % refFlip = (refFlipInt /= 0)
      else
        self % refFlip = .false.
      end if
      if (self % diagEnabled .and. (.not. diagFileOpen)) then
        open(unit=DIAG_UNIT, file=trim(diagFile), status='replace', action='write')
        write(DIAG_UNIT, '(A)') '# x  y  z  neural_hs  ref_hs  sdf'
        diagFileOpen = .true.
      end if
    end if

    !! !!TO BE REMOVED!! Optional diagnostic: UNION-of-3 reference surface for
    !! complex multi-shape geometry (the teapot) -- same convention as
    !! deepLSSurface_class.f90's identical addition. Requires all three
    !! of body/handle/spout. Takes precedence over refSurface if both present.
    if (dict % isPresent('refSurfaces')) then
      if (dict % isPresent('diagEnabled')) then
        call dict % get(diagEnabledInt, 'diagEnabled')
        self % diagEnabled = (diagEnabledInt /= 0)
      else
        self % diagEnabled = .true.
      end if
      if (associated(self % refSurfBody))   call self % refSurfBody   % kill()
      if (associated(self % refSurfHandle)) call self % refSurfHandle % kill()
      if (associated(self % refSurfSpout))  call self % refSurfSpout  % kill()

      call dict % get(refSurfacesDict, 'refSurfaces')
      if (.not. refSurfacesDict % isPresent('body'))   call fatalError(Here, "refSurfaces missing 'body'")
      if (.not. refSurfacesDict % isPresent('handle')) call fatalError(Here, "refSurfaces missing 'handle'")
      if (.not. refSurfacesDict % isPresent('spout'))  call fatalError(Here, "refSurfaces missing 'spout'")

      call refSurfacesDict % get(refSurfDict, 'body')
      allocate(bezierVolume :: self % refSurfBody)
      call self % refSurfBody % init(refSurfDict)

      call refSurfacesDict % get(refSurfDict, 'handle')
      allocate(bezierVolume :: self % refSurfHandle)
      call self % refSurfHandle % init(refSurfDict)

      call refSurfacesDict % get(refSurfDict, 'spout')
      allocate(bezierVolume :: self % refSurfSpout)
      call self % refSurfSpout % init(refSurfDict)

      if (dict % isPresent('refFlip')) then
        call dict % get(refFlipInt, 'refFlip')
        self % refFlip = (refFlipInt /= 0)
      else
        self % refFlip = .false.
      end if
      if (self % diagEnabled .and. (.not. diagFileOpen)) then
        open(unit=DIAG_UNIT, file=trim(diagFile), status='replace', action='write')
        write(DIAG_UNIT, '(A)') '# x  y  z  neural_hs  ref_hs  sdf'
        diagFileOpen = .true.
      end if
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
  !! misclassification diagnostics when a reference surface is configured.
  !!
  !! See surface_inter for details
  !!
  function halfspace(self, r, u) result(hs)
    class(neuralSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    logical(defBool)                        :: ref_hs  !! !!TO BE REMOVED!!
    real(defReal)                           :: c
    c = self % mlp % evaluate(r / self % geomScale) * self % geomScale
    if (abs(c) < self % surfTol()) then
      hs = self % going(r, u)
    else
      hs = c > ZERO
    end if

    !! !!TO BE REMOVED!! Misclassification diagnostic: compare with analytic sphere
    if (self % diagEnabled .and. self % diagRefR2 > ZERO) then
      ref_hs = (r(1)*r(1) + r(2)*r(2) + r(3)*r(3) - self % diagRefR2) > ZERO
      !$omp atomic
      neuralSurf_nCalls = neuralSurf_nCalls + 1_longInt
      if (hs .neqv. ref_hs) then
        !$omp atomic
        neuralSurf_nMisclass = neuralSurf_nMisclass + 1_longInt
        if (diagFileOpen) then
          !$omp critical(neuralDiag)
          write(DIAG_UNIT, '(3ES16.8, 2L3, ES16.8)') r(1), r(2), r(3), hs, ref_hs, c
          !$omp end critical(neuralDiag)
        end if
      end if
    end if

    !! !!TO BE REMOVED!! Misclassification diagnostic: union-of-3 reference
    !! (teapot) takes precedence over the single-surface reference if both
    !! are somehow populated -- same convention as
    !! deepLSSurface_class.f90's identical addition.
    if (self % diagEnabled .and. associated(self % refSurfBody)) then
      ref_hs = self % refSurfBody   % halfspace(r, u) .and. &
               self % refSurfHandle % halfspace(r, u) .and. &
               self % refSurfSpout  % halfspace(r, u)
      if (self % refFlip) ref_hs = .not. ref_hs
      !$omp atomic
      neuralSurf_nCalls = neuralSurf_nCalls + 1_longInt
      if (hs .neqv. ref_hs) then
        !$omp atomic
        neuralSurf_nMisclass = neuralSurf_nMisclass + 1_longInt
        if (diagFileOpen) then
          !$omp critical(neuralDiag)
          write(DIAG_UNIT, '(3ES16.8, 2L3, ES16.8)') r(1), r(2), r(3), hs, ref_hs, c
          !$omp end critical(neuralDiag)
        end if
      end if
    !! !!TO BE REMOVED!! Misclassification diagnostic: compare with reference surface
    else if (self % diagEnabled .and. associated(self % refSurf)) then
      ref_hs = self % refSurf % halfspace(r, u)
      if (self % refFlip) ref_hs = .not. ref_hs
      !$omp atomic
      neuralSurf_nCalls = neuralSurf_nCalls + 1_longInt
      if (hs .neqv. ref_hs) then
        !$omp atomic
        neuralSurf_nMisclass = neuralSurf_nMisclass + 1_longInt
        if (diagFileOpen) then
          !$omp critical(neuralDiag)
          write(DIAG_UNIT, '(3ES16.8, 2L3, ES16.8)') r(1), r(2), r(3), hs, ref_hs, c
          !$omp end critical(neuralDiag)
        end if
      end if
    end if

  end function halfspace

  !!
  !! Print misclassification diagnostic summary  !! !!TO BE REMOVED!!
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
    !! !!TO BE REMOVED!! Release diagnostic reference surface(s)
    if (associated(self % refSurf)) then
      call self % refSurf % kill()
      deallocate(self % refSurf)
    end if
    if (associated(self % refSurfBody)) then
      call self % refSurfBody % kill()
      deallocate(self % refSurfBody)
    end if
    if (associated(self % refSurfHandle)) then
      call self % refSurfHandle % kill()
      deallocate(self % refSurfHandle)
    end if
    if (associated(self % refSurfSpout)) then
      call self % refSurfSpout % kill()
      deallocate(self % refSurfSpout)
    end if

  end subroutine finaliseNeuralSurface

end module neuralSurface_class
