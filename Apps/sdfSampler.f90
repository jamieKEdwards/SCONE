!!
!! SDF training data sampler
!!
!! Samples a SCONE surface to produce (x, y, z, sdf) training data for NeuralSDF.
!! Writes a binary file readable by scripts/neuralSurface/sampler.py:load_scone_binary().
!!
!! Mixed sampling strategy (DeepLS Sec. 4.3):
!!   nearFraction  of nSamples — uniform in bbox, keep only |sdf| < nearDelta
!!   remainder     of nSamples — uniform in bbox (volumetric)
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
!!   nSamples     100000;
!!   nearFraction 0.6;
!!   nearDelta    0.5;
!!   seed         42;
!!   outputFile   sphere_train.bin;
!!   bboxMin      (-7.5 -7.5 -7.5);
!!   bboxMax      ( 7.5  7.5  7.5);
!!
!!   surface {
!!     type sphere;
!!     id 1;
!!     origin (0.0 0.0 0.0);
!!     radius 5.0;
!!   }
!!
!! Notes:
!!   - nearDelta defaults to 0.05 * diagonal of bbox if not specified.
!!   - seed defaults to 1 if not specified.
!!   - nearFraction defaults to 0.5 if not specified.
!!   - Records contain surface%evaluate(r), which is the implicit equation
!!     residual F(r), NOT a proper signed distance in physical units.
!!     For SCONE's quadratic surfaces (sphere, cylinder, ...) F(r) shares
!!     the same zero set and sign as the true SDF, which is all that delta
!!     (Woodcock) tracking requires. The trained MLP reproduces F(r).
!!   - Use Python train.py --input <outputFile> to train.
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
  type(dictionary)       :: input
  type(dictionary)       :: surfDict
  class(surface), pointer :: surf => null()
  type(rng)              :: rand
  character(:), allocatable :: inputPath
  character(pathLen)        :: outputFile
  integer(shortInt)  :: nSamples, seed_si
  real(defReal)      :: nearFraction, nearDelta
  real(defReal), dimension(3) :: bboxMin, bboxMax
  real(defReal), dimension(:), allocatable :: bbMin3, bbMax3
  character(100), parameter :: Here = 'sdfSampler (sdfSampler.f90)'

  ! ---- Sampling buffers ----
  real(defReal), allocatable :: volBuf(:,:)   ! (4, nVol)   — volumetric points
  real(defReal), allocatable :: nearBuf(:,:)  ! (4, nNear)  — near-surface points
  integer(shortInt) :: nVol, nNear
  integer(shortInt) :: nVolDone, nNearDone
  integer(shortInt) :: batchSize, nBatch
  real(defReal)     :: x, y, z, sdfVal
  real(defReal), dimension(3) :: r
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
    ! Default: 5% of bbox diagonal
    nearDelta = 0.05_defReal * norm2(bboxMax - bboxMin)
  end if

  ! ---- Build surface ----
  call input % get(surfDict, 'surface')
  surf => new_surface_ptr(surfDict)

  ! ---- Initialise RNG ----
  call rand % init(int(seed_si, int64))

  ! ---- Determine sample counts ----
  nNear = int(nSamples * nearFraction, shortInt)
  nVol  = nSamples - nNear

  allocate(volBuf (4, nVol ))
  allocate(nearBuf(4, nNear))

  ! ---- Volumetric sampling ----
  do i = 1, nVol
    r(1) = bboxMin(1) + rand % get() * (bboxMax(1) - bboxMin(1))
    r(2) = bboxMin(2) + rand % get() * (bboxMax(2) - bboxMin(2))
    r(3) = bboxMin(3) + rand % get() * (bboxMax(3) - bboxMin(3))
    volBuf(1, i) = r(1)
    volBuf(2, i) = r(2)
    volBuf(3, i) = r(3)
    volBuf(4, i) = surf % evaluate(r)
  end do

  ! ---- Near-surface sampling (rejection) ----
  ! Oversample in batches of batchSize until nNear accepted samples collected
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

  ! ---- Write binary output ----
  call writeBinary(trim(outputFile), nVol, nNear, volBuf, nearBuf, bboxMin, bboxMax)

  ! ---- Clean up ----
  call surf % kill()
  deallocate(surf)
  deallocate(volBuf)
  deallocate(nearBuf)
  call input % kill()
  call surfDict % kill()

  print '(A,I0,A)', 'sdfSampler: wrote ', nVol + nNear, ' records to '//trim(outputFile)
  print '(A,I0,A,I0)', '  volumetric: ', nVol, '   near-surface: ', nNear
  print '(A,F8.4)', '  nearDelta: ', nearDelta

contains

  !!
  !! Write the SCONE sdfSampler binary format.
  !!
  !! Format:
  !!   count   int32
  !!   bbox    float64 × 6   (xmin ymin zmin xmax ymax zmax)
  !!   records float64 × 4   repeated count times (x y z sdf)
  !!
  subroutine writeBinary(filename, nVol, nNear, volBuf, nearBuf, bboxMin, bboxMax)
    character(*), intent(in)              :: filename
    integer(shortInt), intent(in)         :: nVol, nNear
    real(defReal), intent(in)             :: volBuf(4, nVol)
    real(defReal), intent(in)             :: nearBuf(4, nNear)
    real(defReal), dimension(3), intent(in) :: bboxMin, bboxMax
    integer(shortInt)  :: i, unit
    integer(int32)     :: count32
    real(8)            :: bbox8(6), rec8(4)
    character(100), parameter :: Here = 'writeBinary (sdfSampler.f90)'

    count32 = int(nVol + nNear, int32)
    bbox8(1:3) = real(bboxMin, 8)
    bbox8(4:6) = real(bboxMax, 8)

    open(newunit=unit, file=filename, status='replace', &
         access='stream', form='unformatted', iostat=i)
    if (i /= 0) call fatalError(Here, 'Cannot open output file: '//trim(filename))

    ! Header
    write(unit) count32
    write(unit) bbox8

    ! Volumetric records
    do i = 1, nVol
      rec8 = real(volBuf(:, i), 8)
      write(unit) rec8
    end do

    ! Near-surface records
    do i = 1, nNear
      rec8 = real(nearBuf(:, i), 8)
      write(unit) rec8
    end do

    close(unit)

  end subroutine writeBinary

end program sdfSampler
