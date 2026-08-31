module deepLSSurface_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill
  use trainedMLP_class,   only : trainedMLP, normaliseToBox

  implicit none
  private

  character(*), parameter :: TYPE_NAME = 'deepLSSurface'

  !! Voxel status codes (must match export_deepls.py's STATUS_* constants)
  integer(shortInt), parameter :: DEEPLS_CONST_INSIDE  = 0_shortInt
  integer(shortInt), parameter :: DEEPLS_CONST_OUTSIDE = 1_shortInt
  integer(shortInt), parameter :: DEEPLS_HAS_MLP       = 2_shortInt

  !! Weight file format (see scripts/neuralSurface/export_deepls.py). The magic
  !! bytes are the same "NSDF" tag as the single-MLP file (trainedMLP_class),
  !! but the layout past the header is different: version 3, with a voxel grid.
  integer(shortInt), parameter :: MAGIC_NUMBER   = int(z'4E534446', shortInt)
  integer(shortInt), parameter :: FORMAT_VERSION = 3_shortInt
  real(defReal),     parameter :: VALIDATION_TOL = 1.0e-9_defReal

  !!
  !! DeepLS surface, paper-faithful architecture (Chabra et al., ECCV 2020)
  !!
  !! ONE shared decoder network (same size as the global MLP: 4 layers, 128
  !! hidden by default) conditioned on a per-voxel latent code, rather than one
  !! global network representing the whole domain. The decoder is evaluated
  !! identically for every voxel; only the latent code and the local
  !! normalisation bbox change. Trained via scripts/neuralSurface/train_deepls.py.
  !!
  !! Sign convention identical to neuralSurface: F(r) < 0 => inside,
  !! F(r) > 0 => outside. Delta (Woodcock) tracking only.
  !!
  !! See misclassClerk_class (Tallies/TallyClerks) for a halfspace
  !! misclassification diagnostic against a reference region.
  !!
  !! Private Members:
  !!   decoder      -> the ONE shared trainedMLP (inputDim = latentDim + 3)
  !!   latentDim    -> latent code length
  !!   nvox         -> voxel grid dimensions
  !!   gridOrigin   -> world-space min corner of the voxel grid
  !!   voxelSize    -> world-space size of one voxel, per axis
  !!   voxelStatus  -> (nx,ny,nz) DEEPLS_* status codes
  !!   voxelBboxMin -> (3,nx,ny,nz) per-voxel local normalisation box, min corner
  !!   voxelBboxMax -> (3,nx,ny,nz) per-voxel local normalisation box, max corner
  !!   voxelLatent  -> (latentDim,nx,ny,nz) per-voxel latent code
  !!   geomScale    -> divides world-space coordinates before the voxel-grid
  !!                   lookup, letting one trained weight file represent
  !!                   geometrically similar surfaces of a different physical
  !!                   size without retraining
  !!
  !! Interface:
  !!   myType      -> Return surface type name
  !!   init        -> Build from a dictionary (reads the weight file)
  !!   boundingBox -> Voxel-grid extent scaled by geomScale
  !!   evaluate    -> Signed surface value at a point
  !!   distance    -> Stub (INF); delta-tracking-only surface
  !!   going       -> Finite-difference halfspace resolution on the boundary
  !!   kill        -> Free the grid arrays and reset to the uninitialised state
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
  contains
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
    procedure, private :: killGrids
    procedure, private :: loadWeights
    procedure, private :: decoderInput
    procedure, private :: validateTestVector
  end type deepLSSurface

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for details
  !!
  pure function myType(self) result(str)
    class(deepLSSurface), intent(in) :: self
    character(:), allocatable        :: str

    str = TYPE_NAME

  end function myType

  !!
  !! Initialise from a dictionary
  !!
  !! Reads 'id' and 'weightFile', loads the DeepLS model from that weight file,
  !! and reads the optional 'geometricScale'.
  !!
  !! See surface_inter for details
  !!
  !! Errors:
  !!   fatalError if id < 1, geometricScale <= 0, or the weight file cannot be read
  !!
  subroutine init(self, dict)
    class(deepLSSurface), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: id
    character(pathLen)                  :: weightFile
    character(100), parameter :: Here = 'init (deepLSSurface_class.f90)'

    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Invalid surface id: '//numToChar(id))

    call dict % get(weightFile, 'weightFile')
    call self % loadWeights(trim(weightFile))
    call self % setId(id)

    ! Optional geometric scale: divides world coordinates before the voxel-grid
    ! lookup, so one weight file can represent similar surfaces of a different
    ! physical size.
    call dict % getOrDefault(self % geomScale, 'geometricScale', ONE)
    if (self % geomScale <= ZERO) call fatalError(Here, 'geometricScale must be positive')

  end subroutine init

  !!
  !! Return the axis-aligned bounding box
  !!
  !! The voxel-grid extent scaled by geomScale.
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(deepLSSurface), intent(in) :: self
    real(defReal), dimension(6)      :: aabb

    aabb(1:3) = self % gridOrigin * self % geomScale
    aabb(4:6) = (self % gridOrigin + self % nvox * self % voxelSize) * self % geomScale

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! Divides r by geomScale, locates the voxel, and returns:
  !!   -1                          for a CONST_INSIDE voxel
  !!   +1                          for a CONST_OUTSIDE voxel or a point outside
  !!                               the grid
  !!   the shared decoder output   for a HAS_MLP voxel
  !! all then multiplied back by geomScale.
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(deepLSSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: rs
    integer(shortInt), dimension(3)         :: idx

    rs = r / self % geomScale

    if (any(rs < self % gridOrigin) .or. &
        any(rs > self % gridOrigin + self % nvox * self % voxelSize)) then
      c = ONE * self % geomScale
      return
    end if

    ! Clamp is load-bearing: a point exactly on the far grid face passes the
    ! check above, and floor() then yields nvox+1.
    idx = floor((rs - self % gridOrigin) / self % voxelSize) + 1
    idx = max(1, min(self % nvox, idx))

    select case (self % voxelStatus(idx(1), idx(2), idx(3)))
      case (DEEPLS_CONST_INSIDE)
        c = -ONE
      case (DEEPLS_CONST_OUTSIDE)
        c = ONE
      case (DEEPLS_HAS_MLP)
        c = self % decoder % evaluateRaw(self % decoderInput(rs, idx(1), idx(2), idx(3)))
      case default
        c = ONE
    end select

    c = c * self % geomScale

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! Stub: always INF. deepLSSurface is for delta (Woodcock) tracking only.
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(deepLSSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d

    d = INF

  end function distance

  !!
  !! Return TRUE if the particle is going into the +ve halfspace
  !!
  !! Uses a single forward finite-difference step along u.
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(hs)
    class(deepLSSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    real(defReal), parameter :: FD_STEP = 1.0e-7_defReal

    hs = self % evaluate(r + FD_STEP * u) > ZERO

  end function going

  !!
  !! Return to the uninitialised state and free the grid arrays
  !!
  !! See surface_inter for details
  !!
  elemental subroutine kill(self)
    class(deepLSSurface), intent(inout) :: self

    call kill_super(self)
    call self % killGrids()

    self % latentDim  = 0
    self % nvox       = 0
    self % gridOrigin = ZERO
    self % voxelSize  = ONE
    self % geomScale  = ONE

  end subroutine kill

  !!
  !! Free the decoder and the per-voxel grid arrays
  !!
  !! Split out from kill() because loadWeights() needs the same release (without
  !! the surface-level reset) when re-loading into an existing surface.
  !!
  elemental subroutine killGrids(self)
    class(deepLSSurface), intent(inout) :: self

    call self % decoder % kill()
    if (allocated(self % voxelStatus))  deallocate(self % voxelStatus)
    if (allocated(self % voxelBboxMin)) deallocate(self % voxelBboxMin)
    if (allocated(self % voxelBboxMax)) deallocate(self % voxelBboxMax)
    if (allocated(self % voxelLatent))  deallocate(self % voxelLatent)

  end subroutine killGrids

  !!
  !! Assemble the shared decoder's input vector for a point in voxel (ix,iy,iz)
  !!
  !! Layout: [ per-voxel latent code (latentDim) , local xyz normalised onto
  !! [-1,1] using that voxel's own bbox (3) ]. Used by both evaluate() at
  !! runtime and validateTestVector()'s self-consistency check.
  !!
  !! Args:
  !!   r        [in] -> point in (geomScale-divided) grid coordinates
  !!   ix/iy/iz [in] -> voxel indices (must be a DEEPLS_HAS_MLP voxel)
  !!
  !! Result:
  !!   Decoder input vector of length latentDim + 3.
  !!
  pure function decoderInput(self, r, ix, iy, iz) result(inp)
    class(deepLSSurface), intent(in)              :: self
    real(defReal), dimension(3), intent(in)       :: r
    integer(shortInt), intent(in)                 :: ix, iy, iz
    real(defReal), dimension(self % latentDim + 3) :: inp

    inp(1 : self % latentDim) = self % voxelLatent(:, ix, iy, iz)
    inp(self % latentDim + 1 : self % latentDim + 3) = &
         normaliseToBox(r, self % voxelBboxMin(:, ix, iy, iz), self % voxelBboxMax(:, ix, iy, iz))

  end function decoderInput

  ! ---------------------------------------------------------------------------
  ! Weight file I/O
  ! ---------------------------------------------------------------------------

  !!
  !! Read a paper-faithful DeepLS model (shared decoder + per-voxel latent
  !! codes) from a binary weight file (v3 format) into this surface
  !!
  !! See scripts/neuralSurface/export_deepls.py's module docstring for the exact
  !! binary layout. There is ONE shared trainedMLP decoder for the whole grid
  !! (inputDim = latentDim + 3); each active voxel carries a small latent-code
  !! vector plus its own local normalisation bbox.
  !!
  !! Args:
  !!   filename [in] -> path to the weight file
  !!
  !! Errors:
  !!   fatalError if the file cannot be opened, has the wrong magic or version,
  !!   is truncated, has an out-of-range active voxel, or fails the embedded
  !!   test-vector check.
  !!
  subroutine loadWeights(self, filename)
    class(deepLSSurface), intent(inout) :: self
    character(*), intent(in)            :: filename
    integer(shortInt)                   :: unit, stat
    integer(shortInt)                   :: magic, version
    integer(shortInt)                   :: hiddenDim, numLayers, activationType
    real(defReal)                       :: leakyAlpha
    integer(shortInt)                   :: nActive
    integer(shortInt)                   :: a, ix, iy, iz, l, inDim, outDim
    integer(shortInt), dimension(3)     :: voxIdx
    real(defReal), dimension(3)         :: bboxMin, bboxMax
    real(defReal), allocatable          :: latentVec(:)
    real(defReal), dimension(3)         :: testInput
    real(defReal)                       :: testOutput
    character(100), parameter :: Here = 'loadWeights (deepLSSurface_class.f90)'

    ! Release anything from a previous load (re-init is allowed)
    call self % killGrids()

    open(newunit = unit, file = filename, access = 'stream', form = 'unformatted', &
         status = 'old', action = 'read', iostat = stat)
    if (stat /= 0) call fatalError(Here, 'Cannot open weight file: '//trim(filename))

    read(unit, iostat=stat) magic;   call ioCheck(stat, 'magic', Here)
    if (magic /= MAGIC_NUMBER) call fatalError(Here, 'Not an NSDF weight file: '//trim(filename))
    read(unit, iostat=stat) version; call ioCheck(stat, 'version', Here)
    if (version /= FORMAT_VERSION) call fatalError(Here, &
      'Expected DeepLS format version '//numToChar(FORMAT_VERSION)// &
      ', got '//numToChar(version)//'.')

    read(unit, iostat=stat) self % latentDim;  call ioCheck(stat, 'latentDim', Here)
    read(unit, iostat=stat) hiddenDim;         call ioCheck(stat, 'hiddenDim', Here)
    read(unit, iostat=stat) numLayers;         call ioCheck(stat, 'numLayers', Here)
    read(unit, iostat=stat) activationType;    call ioCheck(stat, 'activationType', Here)
    read(unit, iostat=stat) leakyAlpha;        call ioCheck(stat, 'leakyAlpha', Here)
    read(unit, iostat=stat) self % nvox(1);    call ioCheck(stat, 'nvox_x', Here)
    read(unit, iostat=stat) self % nvox(2);    call ioCheck(stat, 'nvox_y', Here)
    read(unit, iostat=stat) self % nvox(3);    call ioCheck(stat, 'nvox_z', Here)
    read(unit, iostat=stat) self % gridOrigin; call ioCheck(stat, 'gridOrigin', Here)
    read(unit, iostat=stat) self % voxelSize;  call ioCheck(stat, 'voxelSize', Here)
    read(unit, iostat=stat) nActive;           call ioCheck(stat, 'nActive', Here)

    if (any(self % nvox <= 0)) call fatalError(Here, 'Invalid voxel grid dimensions in file')
    if (self % latentDim <= 0) call fatalError(Here, 'Invalid latentDim in file')

    allocate(self % voxelStatus(self % nvox(1), self % nvox(2), self % nvox(3)))
    read(unit, iostat=stat) self % voxelStatus
    call ioCheck(stat, 'voxelStatus', Here)

    ! Shared decoder: one trainedMLP, inputDim = latentDim + 3
    call self % decoder % init(self % latentDim + 3_shortInt, hiddenDim, numLayers, &
                               activationType, leakyAlpha, ONE, &
                               [ZERO, ZERO, ZERO], [ONE, ONE, ONE])
    do l = 1_shortInt, numLayers
      call self % decoder % layerDims(l, inDim, outDim)
      read(unit, iostat=stat) self % decoder % weights(1:outDim, 1:inDim, l)
      call ioCheck(stat, 'decoder weights', Here)
      read(unit, iostat=stat) self % decoder % biases(1:outDim, l)
      call ioCheck(stat, 'decoder biases', Here)
    end do

    allocate(self % voxelBboxMin(3, self % nvox(1), self % nvox(2), self % nvox(3)))
    allocate(self % voxelBboxMax(3, self % nvox(1), self % nvox(2), self % nvox(3)))
    allocate(self % voxelLatent(self % latentDim, self % nvox(1), self % nvox(2), self % nvox(3)))
    self % voxelBboxMin = ZERO; self % voxelBboxMax = ONE; self % voxelLatent = ZERO
    allocate(latentVec(self % latentDim))

    do a = 1_shortInt, nActive
      read(unit, iostat=stat) voxIdx;    call ioCheck(stat, 'voxel_index', Here)
      read(unit, iostat=stat) bboxMin;   call ioCheck(stat, 'voxel_bboxMin', Here)
      read(unit, iostat=stat) bboxMax;   call ioCheck(stat, 'voxel_bboxMax', Here)
      read(unit, iostat=stat) latentVec; call ioCheck(stat, 'voxel_latent', Here)

      ix = voxIdx(1); iy = voxIdx(2); iz = voxIdx(3)
      if (ix < 1 .or. ix > self % nvox(1) .or. iy < 1 .or. iy > self % nvox(2) .or. &
          iz < 1 .or. iz > self % nvox(3)) &
        call fatalError(Here, 'Active voxel index out of grid bounds')
      if (self % voxelStatus(ix, iy, iz) /= DEEPLS_HAS_MLP) call fatalError(Here, &
        'Active voxel latent block does not match voxel map status')

      self % voxelBboxMin(:, ix, iy, iz) = bboxMin
      self % voxelBboxMax(:, ix, iy, iz) = bboxMax
      self % voxelLatent(:, ix, iy, iz)  = latentVec
    end do

    if (nActive > 0) then
      read(unit, iostat=stat) testInput;  call ioCheck(stat, 'testInput', Here)
      read(unit, iostat=stat) testOutput; call ioCheck(stat, 'testOutput', Here)
      call self % validateTestVector(testInput, testOutput, Here)
    end if

    close(unit)

  end subroutine loadWeights

  !!
  !! Run the shared decoder once at the file's embedded test point and check it
  !! matches the embedded output
  !!
  !! The test point lies in the first active (HAS_MLP) voxel in Fortran
  !! column-major scan order, matching export_deepls.py.
  !!
  subroutine validateTestVector(self, testInput, testOutput, Here)
    class(deepLSSurface), intent(in)        :: self
    real(defReal), dimension(3), intent(in) :: testInput
    real(defReal), intent(in)               :: testOutput
    character(*), intent(in)                :: Here
    integer(shortInt) :: jx, jy, jz
    logical(defBool)  :: found
    real(defReal)     :: computed, err
    character(30)     :: errStr

    found = .false.
    scan: do jz = 1_shortInt, self % nvox(3)
      do jy = 1_shortInt, self % nvox(2)
        do jx = 1_shortInt, self % nvox(1)
          if (self % voxelStatus(jx, jy, jz) == DEEPLS_HAS_MLP) then
            computed = self % decoder % evaluateRaw(self % decoderInput(testInput, jx, jy, jz))
            found = .true.
            exit scan
          end if
        end do
      end do
    end do scan

    if (.not. found) return

    err = abs(computed - testOutput)
    if (err > VALIDATION_TOL) then
      write(errStr, '(ES12.4)') err
      call fatalError(Here, 'DeepLS weight file validation FAILED: embedded test &
                            &vector mismatch. |computed - expected| = '// &
                            trim(adjustl(errStr))//'. Likely cause: byte order, transpose, '// &
                            'or voxel-ordering error.')
    end if

  end subroutine validateTestVector

  !!
  !! Abort with fatalError if a stream I/O status is non-zero, naming the field
  !!
  subroutine ioCheck(stat, field, Here)
    integer(shortInt), intent(in) :: stat
    character(*), intent(in)      :: field, Here

    if (stat /= 0) call fatalError(Here, 'I/O error reading field: '//trim(field))

  end subroutine ioCheck

end module deepLSSurface_class
