module deepLSSurface_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface
  use mlpInference_mod,   only : trainedMLP
  use deepLSWeightIO_mod, only : readDeepLSWeights, &
                                  DEEPLS_CONST_INSIDE, DEEPLS_CONST_OUTSIDE, &
                                  DEEPLS_HAS_MLP
  use sphere_class,       only : sphere        !! !!TO BE REMOVED!! (diagnostic refSurface only)
  use bezierVolume_class, only : bezierVolume  !! !!TO BE REMOVED!! (diagnostic refSurfaces union only)

  implicit none
  private

  character(*), parameter :: TYPE_NAME = 'deepLSSurface'

  !! !!TO BE REMOVED!! Module-level diagnostic counters and file handle.
  integer(longInt), public              :: deepLS_nCalls    = 0_longInt
  integer(longInt), public              :: deepLS_nMisclass = 0_longInt
  integer(shortInt), parameter, private :: DIAG_UNIT = 95
  logical(defBool),  save,      private :: diagFileOpen = .false.

  !!
  !! DeepLS surface, paper-faithful architecture (Chabra et al., ECCV 2020)
  !!
  !! ONE shared decoder network (same size as the global MLP: 4 layers, 128
  !! hidden by default) conditioned on a per-voxel latent code, rather than
  !! one global network representing the whole domain. The decoder is
  !! evaluated identically for every voxel; only the latent code and the
  !! local normalisation bbox change. Trained via
  !! scripts/neuralSurface/train_deepls.py.
  !!
  !! Sign convention identical to neuralSurface: F(r) < 0 => inside,
  !! F(r) > 0 => outside. Delta (Woodcock) tracking only.
  !!
  !! Private Members:
  !!   decoder      -> the ONE shared trainedMLP (inputDim = latentDim + 3)
  !!   latentDim    -> latent code length
  !!   nvox/gridOrigin/voxelSize -> voxel grid geometry
  !!   voxelStatus  -> (nx,ny,nz) DEEPLS_* status codes
  !!   voxelBboxMin/Max -> (3,nx,ny,nz) per-voxel local normalisation bbox
  !!   voxelLatent  -> (latentDim,nx,ny,nz) per-voxel latent code
  !!   geomScale    -> divides world-space coordinates before the voxel-grid
  !!                  lookup, letting one trained weight file represent
  !!                  geometrically similar surfaces of a different physical
  !!                  size without retraining
  !!
  !! Sample Dictionary Input:
  !!   surf { type deepLSSurface; id 1; weightFile ./weights.bin;
  !!          # geometricScale 1.0; #
  !!        }
  !!
  type, public, extends(surface) :: deepLSSurface
    private
    type(trainedMLP)                                 :: decoder
    integer(shortInt)                                :: latentDim  = 0
    integer(shortInt), dimension(3)                  :: nvox       = 0
    real(defReal), dimension(3)                      :: gridOrigin = ZERO
    real(defReal), dimension(3)                      :: voxelSize  = ONE
    integer(shortInt), dimension(:,:,:), allocatable :: voxelStatus
    real(defReal), dimension(:,:,:,:), allocatable   :: voxelBboxMin, voxelBboxMax
    real(defReal), dimension(:,:,:,:), allocatable   :: voxelLatent
    real(defReal)                                    :: geomScale  = ONE
    !! !!TO BE REMOVED!! Misclassification diagnostic members.
    !! refSurf: single reference surface (e.g. sphere). refSurfBody/Handle/
    !! Spout: union-of-3 reference for complex multi-shape geometry (e.g.
    !! the teapot -- body+lid, handle, spout, matching union_teapot's exact
    !! unionCell structure), used instead of refSurf when populated.
    !! Coordinate frame: both paths compare directly in world coordinates
    !! after this surface's own geomScale division (see evaluate()) -- for
    !! the teapot, geomScale=1 (no deployment scaling, the training data was
    !! already at full physical scale), so the reference bezierVolume shapes
    !! (defined in the same physical ctrlPts as union_teapot) and this
    !! surface are already in the same coordinate frame with no scale
    !! factor to reconcile.
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
    final     :: finaliseDeepLSSurface
  end type deepLSSurface

  public :: printDeepLSDiagnostics !! !!TO BE REMOVED!!

contains

  pure function myType(self) result(str)
    class(deepLSSurface), intent(in) :: self
    character(:), allocatable        :: str

    str = TYPE_NAME

  end function myType

  !!
  !! Initialise from dictionary
  !!
  !! See surface_inter for details
  !!
  subroutine init(self, dict)
    class(deepLSSurface), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: id
    character(pathLen)                  :: weightFile
    character(pathLen)                  :: diagFile
    integer(shortInt)                   :: refFlipInt
    integer(shortInt)                   :: diagEnabledInt
    type(dictionary)                    :: refSurfDict
    character(nameLen)                  :: refType
    type(dictionary)                    :: refSurfacesDict
    character(100), parameter :: Here = 'init (deepLSSurface_class.f90)'

    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Invalid surface id: '//numToChar(id))

    call dict % get(weightFile, 'weightFile')

    call readDeepLSWeights(self % decoder, self % latentDim, self % nvox, &
                           self % gridOrigin, self % voxelSize, self % voxelStatus, &
                           self % voxelBboxMin, self % voxelBboxMax, self % voxelLatent, &
                           trim(weightFile))
    call self % setId(id)

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
      diagFile = 'deepls_misclass.dat'
    end if

    !! !!TO BE REMOVED!! Optional diagnostic: reference surface (sphere only)
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
        case default
          call fatalError(Here, 'Unsupported refSurface type for deepLSSurface: '//trim(refType))
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
        write(DIAG_UNIT, '(A)') '# x  y  z  deepls_hs  ref_hs  sdf'
        diagFileOpen = .true.
      end if
    end if

    !! !!TO BE REMOVED!! Optional diagnostic: UNION-of-3 reference surface for
    !! complex multi-shape geometry (the teapot). Requires all three of
    !! body/handle/spout -- this isn't a general N-shape mechanism, it's sized
    !! exactly to union_teapot's structure. Reference "inside" = inside body
    !! OR inside handle OR inside spout, matching union_teapot's unionCell
    !! (`surfaces [-2:-3:-4]`) exactly. Takes precedence over refSurf if both
    !! are present (not expected in practice, but not hard-erroring on it --
    !! this is a diagnostic feature, not a production configuration path).
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
        write(DIAG_UNIT, '(A)') '# x  y  z  deepls_hs  ref_hs  sdf'
        diagFileOpen = .true.
      end if
    end if

  end subroutine init

  pure function boundingBox(self) result(aabb)
    class(deepLSSurface), intent(in) :: self
    real(defReal), dimension(6)      :: aabb

    aabb(1:3) = self % gridOrigin * self % geomScale
    aabb(4:6) = (self % gridOrigin + self % nvox * self % voxelSize) * self % geomScale

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(deepLSSurface), intent(in)  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: rs, xn
    integer(shortInt), dimension(3)         :: idx
    real(defReal), dimension(self % latentDim + 3) :: decoderInput

    rs = r / self % geomScale

    if (any(rs < self % gridOrigin) .or. &
        any(rs > self % gridOrigin + self % nvox * self % voxelSize)) then
      c = ONE * self % geomScale
      return
    end if

    idx = floor((rs - self % gridOrigin) / self % voxelSize) + 1
    idx = max(1, min(self % nvox, idx))

    select case (self % voxelStatus(idx(1), idx(2), idx(3)))
      case (DEEPLS_CONST_INSIDE)
        c = -ONE
      case (DEEPLS_CONST_OUTSIDE)
        c = ONE
      case (DEEPLS_HAS_MLP)
        decoderInput(1:self % latentDim) = self % voxelLatent(:, idx(1), idx(2), idx(3))
        xn = TWO * (rs - self % voxelBboxMin(:, idx(1), idx(2), idx(3))) / &
             (self % voxelBboxMax(:, idx(1), idx(2), idx(3)) - &
              self % voxelBboxMin(:, idx(1), idx(2), idx(3))) - ONE
        decoderInput(self % latentDim + 1 : self % latentDim + 3) = xn
        c = self % decoder % evaluateRaw(decoderInput)
      case default
        c = ONE
    end select

    c = c * self % geomScale

  end function evaluate

  pure function distance(self, r, u) result(d)
    class(deepLSSurface), intent(in)  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d

    d = INF

  end function distance

  pure function going(self, r, u) result(hs)
    class(deepLSSurface), intent(in)  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    real(defReal), parameter :: FD_STEP = 1.0e-7_defReal

    hs = self % evaluate(r + FD_STEP * u) > ZERO

  end function going

  !!
  !! Return true if particle is in +ve halfspace — overrides default to add
  !! diagnostics, same core logic (evaluate + surfTol + going tie-break).
  !!
  function halfspace(self, r, u) result(hs)
    class(deepLSSurface), intent(in)  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    logical(defBool)                        :: ref_hs  !! !!TO BE REMOVED!!
    real(defReal)                           :: c

    c = self % evaluate(r)
    if (abs(c) < self % surfTol()) then
      hs = self % going(r, u)
    else
      hs = c > ZERO
    end if

    !! !!TO BE REMOVED!! Misclassification diagnostic: union-of-3 reference
    !! (teapot) takes precedence over the single-surface reference (sphere)
    !! if both are somehow populated.
    if (self % diagEnabled .and. associated(self % refSurfBody)) then
      ! "outside union" = outside body AND outside handle AND outside spout
      ! (De Morgan's of "inside union" = inside body OR inside handle OR
      ! inside spout), matching union_teapot's unionCell semantics exactly.
      ref_hs = self % refSurfBody   % halfspace(r, u) .and. &
               self % refSurfHandle % halfspace(r, u) .and. &
               self % refSurfSpout  % halfspace(r, u)
      if (self % refFlip) ref_hs = .not. ref_hs
      !$omp atomic
      deepLS_nCalls = deepLS_nCalls + 1_longInt
      if (hs .neqv. ref_hs) then
        !$omp atomic
        deepLS_nMisclass = deepLS_nMisclass + 1_longInt
        if (diagFileOpen) then
          !$omp critical(deepLSDiag)
          write(DIAG_UNIT, '(3ES16.8, 2L3, ES16.8)') r(1), r(2), r(3), hs, ref_hs, c
          !$omp end critical(deepLSDiag)
        end if
      end if
    !! !!TO BE REMOVED!! Misclassification diagnostic: single reference (sphere)
    else if (self % diagEnabled .and. associated(self % refSurf)) then
      ref_hs = self % refSurf % halfspace(r, u)
      if (self % refFlip) ref_hs = .not. ref_hs
      !$omp atomic
      deepLS_nCalls = deepLS_nCalls + 1_longInt
      if (hs .neqv. ref_hs) then
        !$omp atomic
        deepLS_nMisclass = deepLS_nMisclass + 1_longInt
        if (diagFileOpen) then
          !$omp critical(deepLSDiag)
          write(DIAG_UNIT, '(3ES16.8, 2L3, ES16.8)') r(1), r(2), r(3), hs, ref_hs, c
          !$omp end critical(deepLSDiag)
        end if
      end if
    end if

  end function halfspace

  subroutine printDeepLSDiagnostics()
    real(defReal) :: pct

    if (deepLS_nCalls == 0_longInt) return

    pct = 100.0_defReal * real(deepLS_nMisclass, defReal) / real(deepLS_nCalls, defReal)

    print '(A)', ''
    print '(A)', '--- DeepLS (Shared Decoder) Surface Halfspace Diagnostic ---'
    print '(A,I0)', '  Total halfspace calls : ', deepLS_nCalls
    print '(A,I0)', '  Misclassified calls   : ', deepLS_nMisclass
    print '(A,F8.4,A)', '  Misclassification rate: ', pct, ' %'
    if (diagFileOpen) then
      write(DIAG_UNIT, '(A)')       '# ---- Summary ----'
      write(DIAG_UNIT, '(A,I0)')    '# Total calls:   ', deepLS_nCalls
      write(DIAG_UNIT, '(A,I0)')    '# Misclassified: ', deepLS_nMisclass
      write(DIAG_UNIT, '(A,F8.4,A)') '# Rate:          ', pct, ' %'
      close(DIAG_UNIT)
      diagFileOpen = .false.
    end if
    print '(A)', '-------------------------------------------------------------'

  end subroutine printDeepLSDiagnostics

  !!
  !! Finaliser -- releases the diagnostic reference surface(s), if allocated
  !!
  subroutine finaliseDeepLSSurface(self)
    type(deepLSSurface), intent(inout) :: self

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

  end subroutine finaliseDeepLSSurface

end module deepLSSurface_class
