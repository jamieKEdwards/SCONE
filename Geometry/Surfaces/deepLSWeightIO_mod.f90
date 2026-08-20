module deepLSWeightIO_mod

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use mlpInference_mod,  only : trainedMLP

  implicit none
  private

  integer(shortInt), parameter :: MAGIC_NUMBER   = int(z'4E534446', shortInt)
  integer(shortInt), parameter :: FORMAT_VERSION = 3_shortInt

  !! Voxel status codes — must match export_deepls.py's STATUS_* constants
  integer(shortInt), parameter, public :: DEEPLS_CONST_INSIDE  = 0_shortInt
  integer(shortInt), parameter, public :: DEEPLS_CONST_OUTSIDE = 1_shortInt
  integer(shortInt), parameter, public :: DEEPLS_HAS_MLP       = 2_shortInt

  real(defReal), parameter :: VALIDATION_TOL = 1.0e-9_defReal

  public :: readDeepLSWeights

contains

  !!
  !! Read a paper-faithful DeepLS model (shared decoder + per-voxel latent
  !! codes) from a binary weight file (v3 format)
  !!
  !! See scripts/neuralSurface/export_deepls.py's module docstring for
  !! the exact binary layout. Unlike the v2 (independent-MLP) format, there is
  !! ONE shared trainedMLP decoder for the whole grid; each active voxel only
  !! carries a small latent-code vector plus its own local normalisation bbox.
  !!
  !! Args:
  !!   decoder      [out] -> the ONE shared trainedMLP (inputDim = latentDim+3)
  !!   latentDim    [out] -> latent code length
  !!   nvox         [out] -> voxel grid dimensions
  !!   gridOrigin   [out] -> world-space min corner of the voxel grid
  !!   voxelSize    [out] -> world-space size of one voxel, per axis
  !!   voxelStatus  [out] -> allocated (nx,ny,nz) array of DEEPLS_* codes
  !!   voxelBboxMin [out] -> allocated (3,nx,ny,nz), valid where status==HAS_MLP
  !!   voxelBboxMax [out] -> allocated (3,nx,ny,nz), valid where status==HAS_MLP
  !!   voxelLatent  [out] -> allocated (latentDim,nx,ny,nz), valid where status==HAS_MLP
  !!   filename     [in]  -> path to weight file
  !!
  subroutine readDeepLSWeights(decoder, latentDim, nvox, gridOrigin, voxelSize, &
                               voxelStatus, voxelBboxMin, voxelBboxMax, voxelLatent, &
                               filename)
    type(trainedMLP), intent(out)                              :: decoder
    integer(shortInt), intent(out)                             :: latentDim
    integer(shortInt), dimension(3), intent(out)                :: nvox
    real(defReal), dimension(3), intent(out)                    :: gridOrigin, voxelSize
    integer(shortInt), dimension(:,:,:), allocatable, intent(out) :: voxelStatus
    real(defReal), dimension(:,:,:,:), allocatable, intent(out)   :: voxelBboxMin, voxelBboxMax
    real(defReal), dimension(:,:,:,:), allocatable, intent(out)   :: voxelLatent
    character(*), intent(in)          :: filename
    integer                           :: unit, stat
    integer(shortInt)                 :: magic, version
    integer(shortInt)                 :: hiddenDim, numLayers, activationType
    real(defReal)                     :: leakyAlpha
    integer(shortInt)                 :: nActive
    integer(shortInt)                 :: a, ix, iy, iz, l, inDim, outDim
    integer(shortInt), dimension(3)   :: voxIdx
    real(defReal), dimension(3)       :: bboxMin, bboxMax
    real(defReal), allocatable        :: latentVec(:)
    real(defReal), dimension(3)       :: testInput
    real(defReal)                     :: testOutput, computed, err, xn
    integer(shortInt)                 :: d
    real(defReal), allocatable        :: decoderInput(:)
    character(30)                     :: errStr
    character(100), parameter :: Here = 'readDeepLSWeights (deepLSWeightIO_mod.f90)'

    open(newunit = unit, file = filename, access = 'stream', form = 'unformatted', &
         status = 'old', action = 'read', iostat = stat)
    if (stat /= 0) call fatalError(Here, 'Cannot open weight file: '//trim(filename))

    read(unit, iostat=stat) magic;   call ioCheck(stat, 'magic', Here)
    if (magic /= MAGIC_NUMBER) call fatalError(Here, 'Not an NSDF weight file: '//trim(filename))
    read(unit, iostat=stat) version; call ioCheck(stat, 'version', Here)
    if (version /= FORMAT_VERSION) call fatalError(Here, &
      'Expected DeepLS format version '//numToChar(FORMAT_VERSION)// &
      ', got '//numToChar(version)//'.')

    read(unit, iostat=stat) latentDim;      call ioCheck(stat, 'latentDim', Here)
    read(unit, iostat=stat) hiddenDim;      call ioCheck(stat, 'hiddenDim', Here)
    read(unit, iostat=stat) numLayers;      call ioCheck(stat, 'numLayers', Here)
    read(unit, iostat=stat) activationType; call ioCheck(stat, 'activationType', Here)
    read(unit, iostat=stat) leakyAlpha;     call ioCheck(stat, 'leakyAlpha', Here)
    read(unit, iostat=stat) nvox(1);        call ioCheck(stat, 'nvox_x', Here)
    read(unit, iostat=stat) nvox(2);        call ioCheck(stat, 'nvox_y', Here)
    read(unit, iostat=stat) nvox(3);        call ioCheck(stat, 'nvox_z', Here)
    read(unit, iostat=stat) gridOrigin;     call ioCheck(stat, 'gridOrigin', Here)
    read(unit, iostat=stat) voxelSize;      call ioCheck(stat, 'voxelSize', Here)
    read(unit, iostat=stat) nActive;        call ioCheck(stat, 'nActive', Here)

    if (any(nvox <= 0)) call fatalError(Here, 'Invalid voxel grid dimensions in file')
    if (latentDim <= 0) call fatalError(Here, 'Invalid latentDim in file')

    allocate(voxelStatus(nvox(1), nvox(2), nvox(3)))
    read(unit, iostat=stat) voxelStatus
    call ioCheck(stat, 'voxelStatus', Here)

    ! Shared decoder: one trainedMLP, inputDim = latentDim + 3
    call decoder % init(latentDim + 3_shortInt, hiddenDim, numLayers, activationType, &
                        leakyAlpha, ONE, [ZERO, ZERO, ZERO], [ONE, ONE, ONE])
    do l = 1_shortInt, numLayers
      call layerDims(l, numLayers, latentDim + 3_shortInt, hiddenDim, inDim, outDim)
      read(unit, iostat=stat) decoder % weights(1:outDim, 1:inDim, l)
      call ioCheck(stat, 'decoder weights', Here)
      read(unit, iostat=stat) decoder % biases(1:outDim, l)
      call ioCheck(stat, 'decoder biases', Here)
    end do

    allocate(voxelBboxMin(3, nvox(1), nvox(2), nvox(3)))
    allocate(voxelBboxMax(3, nvox(1), nvox(2), nvox(3)))
    allocate(voxelLatent(latentDim, nvox(1), nvox(2), nvox(3)))
    voxelBboxMin = ZERO; voxelBboxMax = ONE; voxelLatent = ZERO
    allocate(latentVec(latentDim))

    do a = 1_shortInt, nActive
      read(unit, iostat=stat) voxIdx;  call ioCheck(stat, 'voxel_index', Here)
      read(unit, iostat=stat) bboxMin; call ioCheck(stat, 'voxel_bboxMin', Here)
      read(unit, iostat=stat) bboxMax; call ioCheck(stat, 'voxel_bboxMax', Here)
      read(unit, iostat=stat) latentVec; call ioCheck(stat, 'voxel_latent', Here)

      ix = voxIdx(1); iy = voxIdx(2); iz = voxIdx(3)
      if (ix < 1 .or. ix > nvox(1) .or. iy < 1 .or. iy > nvox(2) .or. iz < 1 .or. iz > nvox(3)) &
        call fatalError(Here, 'Active voxel index out of grid bounds')
      if (voxelStatus(ix, iy, iz) /= DEEPLS_HAS_MLP) call fatalError(Here, &
        'Active voxel latent block does not match voxel map status')

      voxelBboxMin(:, ix, iy, iz) = bboxMin
      voxelBboxMax(:, ix, iy, iz) = bboxMax
      voxelLatent(:, ix, iy, iz)  = latentVec
    end do

    if (nActive > 0) then
      read(unit, iostat=stat) testInput;  call ioCheck(stat, 'testInput', Here)
      read(unit, iostat=stat) testOutput; call ioCheck(stat, 'testOutput', Here)

      firstActiveLoop: block
        integer(shortInt) :: jx, jy, jz
        logical(defBool)  :: found
        found = .false.
        allocate(decoderInput(latentDim + 3_shortInt))
        do jz = 1_shortInt, nvox(3)
          do jy = 1_shortInt, nvox(2)
            do jx = 1_shortInt, nvox(1)
              if (voxelStatus(jx, jy, jz) == DEEPLS_HAS_MLP) then
                do d = 1_shortInt, latentDim
                  decoderInput(d) = voxelLatent(d, jx, jy, jz)
                end do
                do d = 1_shortInt, 3_shortInt
                  xn = TWO * (testInput(d) - voxelBboxMin(d, jx, jy, jz)) / &
                       (voxelBboxMax(d, jx, jy, jz) - voxelBboxMin(d, jx, jy, jz)) - ONE
                  decoderInput(latentDim + d) = xn
                end do
                computed = decoder % evaluateRaw(decoderInput)
                found = .true.
                exit
              end if
            end do
            if (found) exit
          end do
          if (found) exit
        end do
      end block firstActiveLoop

      err = abs(computed - testOutput)
      if (err > VALIDATION_TOL) then
        write(errStr, '(ES12.4)') err
        call fatalError(Here, 'DeepLS weight file validation FAILED: embedded test &
                              &vector mismatch. |computed - expected| = '// &
                              trim(adjustl(errStr))//'. Likely cause: byte order, transpose, '// &
                              'or voxel-ordering error.')
      end if
    end if

    close(unit)

  end subroutine readDeepLSWeights

  !!
  !! Layer input/output dimensions for the shared decoder (inputDim on layer 1
  !! is latentDim+3, not the usual 3 — same convention otherwise)
  !!
  pure subroutine layerDims(l, numLayers, inputDim, hiddenDim, inDim, outDim)
    integer(shortInt), intent(in)  :: l, numLayers, inputDim, hiddenDim
    integer(shortInt), intent(out) :: inDim, outDim

    if (l == 1_shortInt) then
      inDim  = inputDim
      outDim = hiddenDim
    else if (l < numLayers) then
      inDim  = hiddenDim
      outDim = hiddenDim
    else
      inDim  = hiddenDim
      outDim = 1_shortInt
    end if

  end subroutine layerDims

  !!
  !! Fatal-error on a non-zero stream I/O status, naming the field that failed
  !!
  subroutine ioCheck(stat, field, Here)
    integer, intent(in)      :: stat
    character(*), intent(in) :: field, Here
    if (stat /= 0) call fatalError(Here, 'I/O error reading field: '//trim(field))
  end subroutine ioCheck

end module deepLSWeightIO_mod
