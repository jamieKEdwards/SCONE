!!
!! SDF training data sampler
!!
!! Samples a SCONE surface to produce (x, y, z, sdf) training data for NeuralSDF.
!! Writes a binary file readable by scripts/neuralSurface/sampler.py:load_scone_binary().
!!
!! Two sampling modes — auto-detected from the surface type:
!!
!!   SDF mode (default for analytic surfaces):
!!     Uses surface%evaluate(r) to obtain the implicit-equation residual F(r).
!!     F(r) shares the sign and zero set of the true SDF for all SCONE quadratic
!!     surfaces (sphere, cylinder, plane, ...).
!!     Near-surface oversampling is applied: nearFraction of samples are drawn
!!     with |F(r)| < nearDelta to improve boundary accuracy.
!!
!!   Halfspace mode (auto-selected when evaluate always returns 0, e.g. bezierShape):
!!     Uses surface%halfspace(r, u) to obtain inside/outside; writes +1.0 or -1.0.
!!     Near-surface oversampling is skipped (proximity cannot be inferred from
!!     halfspace alone); all nSamples are drawn uniformly from the bbox.
!!     Use --mode binary in train.py with this output.
!!
!! Binary output format:
!!   [Header]
!!     count   int32         total number of records written
!!     bbox    float64 × 6  xmin ymin zmin xmax ymax zmax
!!   [Records, count times]
!!     x y z sdf   float64 × 4
!!
!! Sample input file (SCONE dict format):
!!
!!   nSamples     200000;
!!   nearFraction 0.5;
!!   nearDelta    0.5;
!!   seed         42;
!!   outputFile   trefoil_train.bin;
!!   bboxMin      (-2.0 -2.0 -200.0);
!!   bboxMax      ( 2.0  2.0  200.0);
!!
!!   surface {
!!     type bezierShape;
!!     id 1;
!!     ...
!!   }
!!
!! Notes:
!!   - nearFraction and nearDelta are ignored in halfspace mode.
!!   - seed defaults to 1 if not specified.
!!   - nearFraction defaults to 0.5 if not specified.
!!   - nearDelta defaults to 0.05 * diagonal of bbox if not specified.
!!   - Use Python train.py --input <outputFile> to train.
!!     Add --mode binary --loss bce for halfspace-mode output.
!!
!! Run:
!!   cmake --build Build --target sdfSampler.out
!!   ./Build/sdfSampler.out input.dict
!!
program sdfSampler

  use numPrecision
  use iso_fortran_env,          only : int32, int64
  use genericProcedures,        only : fatalError
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : fileToDict
  use commandLineUI,            only : getInputFile
  use rng_class,                only : rng
  use surface_inter,            only : surface
  use surfaceFactory_func,      only : new_surface_ptr

  implicit none

  ! ---- Locals ----
  type(dictionary)        :: input
  type(dictionary)        :: surfDict
  class(surface), pointer :: surf => null()
  type(rng)               :: rand
  character(:), allocatable :: inputPath
  character(pathLen)        :: outputFile
  integer(shortInt)  :: nSamples, seed_si
  real(defReal)      :: nearFraction, nearDelta
  real(defReal), dimension(3) :: bboxMin, bboxMax, bboxCentre
  real(defReal), dimension(:), allocatable :: bbMin3, bbMax3
  logical(defBool)   :: useHalfspace
  character(100), parameter :: Here = 'sdfSampler (sdfSampler.f90)'

  ! ---- Sampling buffers ----
  real(defReal), allocatable :: volBuf (:,:)   ! (4, nVol)  — volumetric points
  real(defReal), allocatable :: nearBuf(:,:)   ! (4, nNear) — near-surface points
  integer(shortInt) :: nVol, nNear
  integer(shortInt) :: nNearDone
  integer(shortInt) :: batchSize, nBatch
  real(defReal)     :: x, y, z, sdfVal
  real(defReal), dimension(3) :: r, u_dummy
  logical(defBool)  :: hs
  integer(shortInt) :: i

  ! ---- Get input file from command line ----
  call getInputFile(inputPath)
  call fileToDict(input, inputPath)

  ! ---- Read required parameters ----
  call input % get(nSamples,   'nSamples')
  call input % get(outputFile, 'outputFile')
  call input % get(bbMin3,     'bboxMin')
  call input % get(bbMax3,     'bboxMax')

  if (size(bbMin3) /= 3) call fatalError(Here, 'bboxMin must have 3 entries')
  if (size(bbMax3) /= 3) call fatalError(Here, 'bboxMax must have 3 entries')
  bboxMin = bbMin3
  bboxMax = bbMax3

  ! ---- Read optional parameters ----
  if (input % isPresent('seed')) then
    call input % get(seed_si, 'seed')
  else
    seed_si = 1_shortInt
  end if

  if (input % isPresent('nearFraction')) then
    call input % get(nearFraction, 'nearFraction')
  else
    nearFraction = 0.5_defReal
  end if

  if (input % isPresent('nearDelta')) then
    call input % get(nearDelta, 'nearDelta')
  else
    nearDelta = 0.05_defReal * norm2(bboxMax - bboxMin)
  end if

  ! ---- Build surface ----
  call input % get(surfDict, 'surface')
  surf => new_surface_ptr(surfDict)

  ! ---- Auto-detect sampling mode ----
  ! Probe evaluate() at the bbox centre. Surfaces without a meaningful implicit
  ! equation (e.g. bezierShape) return exactly 0 always; analytic surfaces
  ! return a non-zero residual for any point not exactly on the surface.
  bboxCentre = HALF * (bboxMin + bboxMax)
  useHalfspace = (surf % evaluate(bboxCentre) == ZERO)

  if (useHalfspace) then
    print '(A)', 'sdfSampler: evaluate() returns 0 — switching to halfspace mode'
    print '(A)', '  Writing +/-1.0 labels. Use train.py --mode binary --loss bce.'
    nVol  = nSamples
    nNear = 0
  else
    print '(A)', 'sdfSampler: SDF mode (using surface evaluate)'
    nNear = int(nSamples * nearFraction, shortInt)
    nVol  = nSamples - nNear
  end if

  ! ---- Initialise RNG ----
  call rand % init(int(seed_si, int64))

  allocate(volBuf (4, nVol ))
  if (nNear > 0) allocate(nearBuf(4, nNear))

  ! Dummy direction for halfspace queries (only matters exactly on surface)
  u_dummy = [ONE, ZERO, ZERO]

  ! ---- Volumetric sampling ----
  do i = 1, nVol
    r(1) = bboxMin(1) + rand % get() * (bboxMax(1) - bboxMin(1))
    r(2) = bboxMin(2) + rand % get() * (bboxMax(2) - bboxMin(2))
    r(3) = bboxMin(3) + rand % get() * (bboxMax(3) - bboxMin(3))
    volBuf(1, i) = r(1)
    volBuf(2, i) = r(2)
    volBuf(3, i) = r(3)
    if (useHalfspace) then
      hs = surf % halfspace(r, u_dummy)
      volBuf(4, i) = merge(-ONE, ONE, hs)  ! inside (-ve hs) -> -1, outside -> +1
    else
      volBuf(4, i) = surf % evaluate(r)
    end if
  end do

  ! ---- Near-surface sampling (SDF mode only) ----
  if (nNear > 0) then
    batchSize = max(nNear * 4_shortInt, 1000_shortInt)
    nNearDone = 0
    do while (nNearDone < nNear)
      nBatch = min(batchSize, nNear - nNearDone) * 4_shortInt
      do i = 1, nBatch
        r(1) = bboxMin(1) + rand % get() * (bboxMax(1) - bboxMin(1))
        r(2) = bboxMin(2) + rand % get() * (bboxMax(2) - bboxMin(2))
        r(3) = bboxMin(3) + rand % get() * (bboxMax(3) - bboxMin(3))
        sdfVal = surf % evaluate(r)
        if (abs(sdfVal) < nearDelta) then
          nNearDone = nNearDone + 1
          nearBuf(1, nNearDone) = r(1)
          nearBuf(2, nNearDone) = r(2)
          nearBuf(3, nNearDone) = r(3)
          nearBuf(4, nNearDone) = sdfVal
          if (nNearDone == nNear) exit
        end if
      end do
    end do
  end if

  ! ---- Write binary output ----
  if (nNear > 0) then
    call writeBinary(trim(outputFile), nVol, nNear, volBuf, nearBuf, bboxMin, bboxMax)
  else
    call writeBinary(trim(outputFile), nVol, 0, volBuf, volBuf(:, 1:0), bboxMin, bboxMax)
  end if

  ! ---- Clean up ----
  call surf % kill()
  deallocate(surf)
  deallocate(volBuf)
  if (allocated(nearBuf)) deallocate(nearBuf)
  call input % kill()
  call surfDict % kill()

  print '(A,I0,A)', 'sdfSampler: wrote ', nVol + nNear, ' records to '//trim(outputFile)
  if (useHalfspace) then
    print '(A,I0)', '  halfspace mode — uniform samples: ', nVol
  else
    print '(A,I0,A,I0)', '  volumetric: ', nVol, '   near-surface: ', nNear
    print '(A,F8.4)', '  nearDelta: ', nearDelta
  end if

contains

  subroutine writeBinary(filename, nVol, nNear, volBuf, nearBuf, bboxMin, bboxMax)
    character(*), intent(in)                :: filename
    integer(shortInt), intent(in)           :: nVol, nNear
    real(defReal), intent(in)               :: volBuf(4, nVol)
    real(defReal), intent(in)               :: nearBuf(4, nNear)
    real(defReal), dimension(3), intent(in) :: bboxMin, bboxMax
    integer(shortInt) :: i, unit
    integer(int32)    :: count32
    real(8)           :: bbox8(6), rec8(4)
    character(100), parameter :: Here = 'writeBinary (sdfSampler.f90)'

    count32 = int(nVol + nNear, int32)
    bbox8(1:3) = real(bboxMin, 8)
    bbox8(4:6) = real(bboxMax, 8)

    open(newunit=unit, file=filename, status='replace', &
         access='stream', form='unformatted', iostat=i)
    if (i /= 0) call fatalError(Here, 'Cannot open output file: '//trim(filename))

    write(unit) count32
    write(unit) bbox8
    do i = 1, nVol
      rec8 = real(volBuf(:, i), 8)
      write(unit) rec8
    end do
    do i = 1, nNear
      rec8 = real(nearBuf(:, i), 8)
      write(unit) rec8
    end do
    close(unit)

  end subroutine writeBinary

end program sdfSampler
