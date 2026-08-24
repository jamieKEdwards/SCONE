!!
!! SDF training data sampler
!!
!! Samples a SCONE surface (or union of surfaces) to produce (x, y, z, sdf)
!! training data for NeuralSDF. Writes a binary file readable by
!! scripts/neuralSurface/sampler.py:load_scone_binary().
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
!!   Halfspace mode (auto-selected when evaluate always returns 0, e.g. bezierShape,
!!   bezierVolume; always used for a "surfaces" union):
!!     Uses surface%halfspace(r, u) to obtain inside/outside; writes +1.0 or -1.0.
!!     Near-surface oversampling, when available, is done by sampling directly ON
!!     the surface rather than by rejection sampling: for bezierVolume-format
!!     surfaces (dictionaries with "numPatches"/"ctrlPts"), a random patch and
!!     parametric (u,v) location are drawn, the patch is evaluated there via
!!     rational-Bezier De Casteljau, and the resulting point is offset along its
!!     (finite-difference) normal to give two near-surface samples for two
!!     halfspace() calls, with no rejection loop. This deliberately mirrors the
!!     three-part near-surface strategy in
!!     scripts/neuralSurface/sampler.py:generate_sphere_sdf (same fractions and
!!     distance parameters, applied to Bezier patches instead of an analytic
!!     sphere), split into three bands:
!!       - Tight linear:  nearFraction/2 of samples,
!!         offset ~ Uniform(-nearDistanceTight, nearDistanceTight)
!!       - Loose linear:  nearFraction/2 of samples, offset ~ Uniform(-nearDistance, nearDistance)
!!       - Gaussian:      surfaceFraction of samples, offset ~ N(0, 0.3*nearDistance)
!!     with the remaining 1 - nearFraction - surfaceFraction drawn uniformly
!!     over the bbox (volumetric). The labels always come from the FULL
!!     geometry (the single surface, or the union), so overlap between union
!!     members (e.g. the teapot handle re-entering the body) is still labelled
!!     correctly.
!!     This requires "numPatches"/"ctrlPts" on every surface involved (every union
!!     member, or the single surface); if unavailable (e.g. bezierShape, which has
!!     no patch grid), near-surface oversampling is skipped and all samples are
!!     drawn uniformly instead.
!!     Use --mode binary in train.py with this output.
!!
!! Single-surface vs. union input:
!!   Exactly one of "surface" (a single surface definition) or "surfaces" (named
!!   sub-dictionaries of multiple surfaces, combined as a CSG union — inside iff
!!   inside ANY constituent surface) must be given. A union always uses halfspace
!!   mode, since no single combined evaluate() residual exists across independent
!!   surfaces.
!!
!! Binary output format:
!!   [Header]
!!     count   int32         total number of records written
!!     bbox    float64 × 6  xmin ymin zmin xmax ymax zmax
!!   [Records, count times]
!!     x y z sdf   float64 × 4
!!
!! Sample input file, single surface (SCONE dict format):
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
!! Sample input file, union of surfaces (e.g. a multi-part Bezier assembly):
!!
!!   nSamples          400000;
!!   nearFraction      0.65;
!!   surfaceFraction   0.10;
!!   nearDistance      0.05;
!!   nearDistanceTight 0.005;
!!   seed              42;
!!   outputFile   teapot_train.bin;
!!   bboxMin      (-3.6375 -2.4 -0.325);
!!   bboxMax      ( 4.1625  2.4  3.475);
!!
!!   surfaces {
!!     body   { type bezierVolume; id 2; numPatches ...; ctrlPts (...); }
!!     handle { type bezierVolume; id 3; numPatches ...; ctrlPts (...); }
!!     spout  { type bezierVolume; id 4; numPatches ...; ctrlPts (...); }
!!   }
!!
!! Notes:
!!   - nearFraction means something different in the two modes: in SDF mode it
!!     is the total near-surface fraction; in halfspace mode it is the
!!     combined tight+loose linear-band fraction (split evenly between them),
!!     matching sampler.py's near_fraction semantics exactly -- the on-surface
!!     Gaussian band is controlled separately by surfaceFraction.
!!   - nearDelta/nearDistance/nearDistanceTight/surfaceFraction are ignored in
!!     SDF mode; nearDelta is ignored in halfspace mode; surfaceFraction/
!!     nearDistance/nearDistanceTight are ignored in halfspace mode without
!!     bezierVolume-format patch data available.
!!   - seed defaults to 1 if not specified.
!!   - nearFraction defaults to 0.5 in SDF mode, 0.65 in halfspace mode.
!!   - nearDelta defaults to 0.05 * diagonal of bbox if not specified.
!!   - surfaceFraction (halfspace mode only) defaults to 0.10 if not specified.
!!   - nearDistance (halfspace mode only; loose-band bound, and 0.3*nearDistance
!!     is the Gaussian band's standard deviation) defaults to 0.05 * diagonal
!!     of bbox if not specified.
!!   - nearDistanceTight (halfspace mode only; tight-band bound) defaults to
!!     0.1 * nearDistance if not specified.
!!   - Use Python train.py --input <outputFile> to train.
!!     Add --mode binary --loss bce for halfspace-mode output.
!!
!! Run:
!!   cmake --build Build --target sdfSampler.out
!!   ./Build/sdfSampler.out input.dict
!!
program sdfSampler

  use numPrecision
  use iso_fortran_env,          only : int32, int64, real64
  use genericProcedures,        only : fatalError, numToChar, crossProduct
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : fileToDict
  use commandLineUI,            only : getInputFile
  use rng_class,                only : rng
  use surface_inter,            only : surface
  use surfaceFactory_func,      only : new_surface_ptr
  use bezierPatch_func,         only : evalPatch

  implicit none

  !!
  !! Small, local container to store a pointer to a polymorphic surface in an array
  !!
  !! Public Members:
  !!   ptr -> Pointer to a surface
  !!
  type :: surfPtr
    class(surface), pointer :: ptr => null()
  end type surfPtr

  !!
  !! Raw bicubic Bezier patch data for one bezierVolume-format surface, parsed
  !! directly from its dictionary (the surface object's own ctrlPts/weights are
  !! private, so this is duplicated here to allow direct on-surface sampling).
  !!
  !! Public Members:
  !!   ctrlPts    -> (numPatches, 4, 4, 3) control points; ctrlPts(p,i,j,:) has
  !!     i = u-index, j = v-index, matching bezierVolume_class's convention
  !!   weights    -> (numPatches, 4, 4) rational weights (all 1.0 if the surface
  !!     dictionary has no "weights" entry, i.e. polynomial Bezier)
  !!   numPatches -> Number of patches
  !!
  type :: patchData
    real(defReal), dimension(:,:,:,:), allocatable :: ctrlPts
    real(defReal), dimension(:,:,:), allocatable   :: weights
    integer(shortInt)                              :: numPatches = 0
  end type patchData

  !!
  !! Bundle of the surface(s) being sampled, used to keep the near-surface
  !! sampling call signatures below short
  !!
  !! Public Members:
  !!   useUnion -> True if surfs (union) should be used instead of surf
  !!   surf     -> Single-surface pointer, used when useUnion is false
  !!   surfs    -> Array of surfaces making up the union, used when useUnion is true
  !!
  type :: sampleTarget
    logical(defBool)                         :: useUnion
    class(surface), pointer                  :: surf
    type(surfPtr), dimension(:), allocatable :: surfs
  end type sampleTarget

  ! ---- Default values ----

  ! General
  integer(shortInt), parameter :: DEFAULT_SEED          = 1_shortInt
  real(defReal),     parameter :: DEFAULT_NEAR_FRACTION = 0.5_defReal ! SDF mode default

  ! Rejection-sampling batch sizing, used by SDF-mode near-surface sampling
  integer(shortInt), parameter :: NEAR_BATCH_MULTIPLIER = 4_shortInt
  integer(shortInt), parameter :: MIN_BATCH_SIZE         = 1000_shortInt

  ! Retry cap for degenerate-normal / non-converged on-surface sampling attempts
  integer(shortInt), parameter :: DEGENERACY_RETRY_MULTIPLIER = 20_shortInt
  integer(shortInt), parameter :: MIN_DEGENERACY_RETRIES      = 10000_shortInt

  ! Finite-difference step in (u,v) patch parameter space, used to estimate
  ! the surface normal for direct on-surface near-surface sampling
  real(defReal), parameter :: FD_PARAM_EPS = 1.0E-4_defReal

  ! Degeneracy guards: a near-zero surface normal, and a near-zero Box-Muller
  ! log() argument
  real(defReal), parameter :: MIN_NORMAL_LENGTH = 1.0E-8_defReal
  real(defReal), parameter :: MIN_LOG_ARG       = 1.0E-300_defReal

  ! Halfspace-mode near-surface sampling defaults. DEFAULT_NEAR_FRACTION_HS and
  ! DEFAULT_SURFACE_FRACTION match scripts/neuralSurface/sampler.py's near_fraction
  ! and on-surface fractions; nearDelta/nearDistance both default to
  ! DEFAULT_PROXIMITY_FACTOR * bbox diagonal; nearDistanceTight defaults to
  ! DEFAULT_TIGHT_LOOSE_RATIO * nearDistance; the on-surface Gaussian band's
  ! standard deviation is GAUSSIAN_SIGMA_FACTOR * nearDistance, matching
  ! sampler.py's hardcoded on-surface perturbation scale.
  real(defReal), parameter :: DEFAULT_NEAR_FRACTION_HS  = 0.65_defReal
  real(defReal), parameter :: DEFAULT_SURFACE_FRACTION  = 0.10_defReal
  real(defReal), parameter :: DEFAULT_PROXIMITY_FACTOR  = 0.05_defReal
  real(defReal), parameter :: DEFAULT_TIGHT_LOOSE_RATIO = 0.1_defReal
  real(defReal), parameter :: GAUSSIAN_SIGMA_FACTOR     = 0.3_defReal

  ! ---- Locals ----
  type(dictionary)                              :: input
  type(dictionary)                              :: surfDict
  type(dictionary)                              :: surfsDict
  class(dictionary), pointer                    :: dPtr
  class(surface), pointer                       :: surf
  type(surfPtr), dimension(:), allocatable      :: surfs
  type(sampleTarget)                            :: tgt
  type(patchData), dimension(:), allocatable    :: patches
  character(nameLen), dimension(:), allocatable :: surfNames
  logical(defBool)                              :: haveSingle, haveUnion, useUnion
  logical(defBool)                              :: havePatches, haveEval
  type(rng)                                     :: rand
  character(:), allocatable                     :: inputPath
  character(pathLen)                            :: outputFile
  integer(shortInt)                             :: nSamples, seed_si
  real(defReal)                                 :: nearFraction, nearDelta
  real(defReal)                                 :: nearDistance, nearDistanceTight, surfaceFraction
  logical(defBool)                              :: nearFractionGiven
  real(defReal), dimension(3)                   :: bboxMin, bboxMax, bboxCentre
  real(defReal), dimension(:), allocatable      :: bbMin3, bbMax3
  logical(defBool)                              :: useHalfspace
  character(100), parameter :: Here = 'sdfSampler (sdfSampler.f90)'

  ! ---- Sampling buffers ----
  real(defReal), allocatable  :: volBuf (:,:)   ! (4, nVol)  — volumetric points
  real(defReal), allocatable  :: nearBuf(:,:)   ! (4, nNear) — near-surface points
  integer(shortInt)           :: nVol, nNear
  integer(shortInt)           :: nTight, nLoose, nGauss
  integer(shortInt)           :: nNearDone
  integer(shortInt)           :: batchSize, nBatch
  real(defReal)               :: sdfVal
  real(defReal), dimension(3) :: r, uDummy
  logical(defBool)            :: hs
  integer(shortInt)           :: i

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
    seed_si = DEFAULT_SEED
  end if

  ! nearFraction's default differs by mode (SDF vs halfspace); its default is
  ! applied later, once useHalfspace is known.
  nearFractionGiven = input % isPresent('nearFraction')
  if (nearFractionGiven) call input % get(nearFraction, 'nearFraction')

  if (input % isPresent('nearDelta')) then
    call input % get(nearDelta, 'nearDelta')
  else
    nearDelta = DEFAULT_PROXIMITY_FACTOR * norm2(bboxMax - bboxMin)
  end if

  if (input % isPresent('surfaceFraction')) then
    call input % get(surfaceFraction, 'surfaceFraction')
  else
    surfaceFraction = DEFAULT_SURFACE_FRACTION
  end if

  if (input % isPresent('nearDistance')) then
    call input % get(nearDistance, 'nearDistance')
  else
    nearDistance = DEFAULT_PROXIMITY_FACTOR * norm2(bboxMax - bboxMin)
  end if

  if (input % isPresent('nearDistanceTight')) then
    call input % get(nearDistanceTight, 'nearDistanceTight')
  else
    nearDistanceTight = DEFAULT_TIGHT_LOOSE_RATIO * nearDistance
  end if

  ! ---- Build surface(s) ----
  haveSingle = input % isPresent('surface')
  haveUnion  = input % isPresent('surfaces')

  if (haveSingle .eqv. haveUnion) then
    call fatalError(Here, 'Input must contain exactly one of "surface" (single &
                          &surface) or "surfaces" (union of surfaces), not both &
                          &or neither.')
  end if

  useUnion = haveUnion
  allocate(surfs(0))

  if (useUnion) then
    call input % get(surfsDict, 'surfaces')
    call surfsDict % keys(surfNames, 'dict')

    if (size(surfNames) < 1) then
      call fatalError(Here, 'surfaces block must contain at least one surface &
                            &sub-dictionary')
    end if

    deallocate(surfs)
    allocate(surfs(size(surfNames)))
    do i = 1, size(surfNames)
      surfs(i) % ptr => new_surface_ptr(surfsDict % getDictPtr(surfNames(i)))
    end do
  else
    call input % get(surfDict, 'surface')
    surf => new_surface_ptr(surfDict)
  end if

  tgt % useUnion = useUnion
  if (useUnion) then
    tgt % surfs = surfs
  else
    tgt % surf => surf
  end if

  ! ---- Auto-detect sampling mode ----
  ! Probe evaluate() at the bbox centre. Surfaces without a meaningful implicit
  ! equation (e.g. bezierShape, bezierVolume) return exactly 0 always; analytic
  ! surfaces return a non-zero residual for any point not exactly on the surface.
  ! A union always uses halfspace mode.
  bboxCentre = HALF * (bboxMin + bboxMax)
  if (useUnion) then
    useHalfspace = .true.
  else
    useHalfspace = (surf % evaluate(bboxCentre) == ZERO)
  end if

  ! ---- Check whether every surface involved has bezierVolume-format patch data,
  !      needed for direct on-surface near-surface sampling ----
  havePatches = .false.
  if (useHalfspace) then
    if (useUnion) then
      havePatches = .true.
      do i = 1, size(surfNames)
        dPtr => surfsDict % getDictPtr(surfNames(i))
        if (.not. dPtr % isPresent('numPatches')) havePatches = .false.
      end do
    else
      havePatches = surfDict % isPresent('numPatches')
    end if
  end if

  ! ---- If no patch data, check whether every surface has a genuine (non-dummy)
  !      evaluate() instead, enabling generic Newton-projection on-surface
  !      near-surface sampling (sampleEvalPoint) as a fallback to the
  !      patch-based method. A surface without a real implicit equation (e.g.
  !      bezierVolume, bezierShape) always returns exactly 0 regardless of
  !      point, exactly like the top-level mode-detection probe above.
  haveEval = .false.
  if (useHalfspace .and. .not. havePatches) then
    if (useUnion) then
      haveEval = .true.
      do i = 1, size(surfs)
        if (surfs(i) % ptr % evaluate(bboxCentre) == ZERO) haveEval = .false.
      end do
    else
      haveEval = (surf % evaluate(bboxCentre) /= ZERO)
    end if
  end if

  allocate(patches(0))
  if (havePatches) then
    if (useUnion) then
      deallocate(patches)
      allocate(patches(size(surfNames)))
      do i = 1, size(surfNames)
        call readPatchData(surfsDict % getDictPtr(surfNames(i)), patches(i))
      end do
    else
      deallocate(patches)
      allocate(patches(1))
      call readPatchData(surfDict, patches(1))
    end if
  end if

  if (useHalfspace) then
    if (useUnion) then
      print '(A)', 'sdfSampler: halfspace mode (union of '//numToChar(size(surfs))//' surfaces)'
    else
      print '(A)', 'sdfSampler: evaluate() returns 0 — switching to halfspace mode'
    end if
    print '(A)', '  Writing +/-1.0 labels. Use train.py --mode binary --loss bce.'
    if (haveEval) then
      print '(A)', '  No bezierVolume-format patch data, but every surface has a real &
                    &evaluate() -- using generic Newton-projection near-surface sampling.'
    else if (.not. havePatches) then
      print '(A)', '  No bezierVolume-format patch data (numPatches/ctrlPts) and no usable &
                    &evaluate() on every surface -- near-surface oversampling disabled, &
                    &sampling uniformly.'
    end if
  else
    print '(A)', 'sdfSampler: SDF mode (using surface evaluate)'
  end if

  if (.not. nearFractionGiven) then
    if (useHalfspace) then
      nearFraction = DEFAULT_NEAR_FRACTION_HS
    else
      nearFraction = DEFAULT_NEAR_FRACTION
    end if
  end if

  if (useHalfspace .and. (havePatches .or. haveEval)) then
    nTight = int(nSamples * nearFraction * HALF, shortInt)
    nLoose = int(nSamples * nearFraction * HALF, shortInt)
    nGauss = int(nSamples * surfaceFraction, shortInt)
    nNear  = nTight + nLoose + nGauss
    nVol   = nSamples - nNear
  else if (useHalfspace) then
    nTight = 0
    nLoose = 0
    nGauss = 0
    nNear  = 0
    nVol   = nSamples
  else
    nTight = 0
    nLoose = 0
    nGauss = 0
    nNear  = int(nSamples * nearFraction, shortInt)
    nVol   = nSamples - nNear
  end if

  ! ---- Initialise RNG ----
  call rand % init(int(seed_si, int64))

  allocate(volBuf (4, nVol ))
  if (nNear > 0) allocate(nearBuf(4, nNear))

  ! Dummy direction for halfspace queries (only matters exactly on surface)
  uDummy = [ONE, ZERO, ZERO]

  ! ---- Volumetric sampling ----
  do i = 1, nVol
    r(1) = bboxMin(1) + rand % get() * (bboxMax(1) - bboxMin(1))
    r(2) = bboxMin(2) + rand % get() * (bboxMax(2) - bboxMin(2))
    r(3) = bboxMin(3) + rand % get() * (bboxMax(3) - bboxMin(3))
    volBuf(1, i) = r(1)
    volBuf(2, i) = r(2)
    volBuf(3, i) = r(3)
    if (useHalfspace) then
      if (useUnion) then
        hs = unionHalfspace(surfs, r, uDummy)
      else
        hs = surf % halfspace(r, uDummy)
      end if
      volBuf(4, i) = merge(-ONE, ONE, hs)  ! inside (-ve hs) -> -1, outside -> +1
    else
      volBuf(4, i) = surf % evaluate(r)
    end if
  end do

  ! ---- Near-surface sampling ----
  if (nNear > 0) then
    if (useHalfspace .and. havePatches) then
      call sampleNearSurf(nTight, nLoose, nGauss, tgt, patches, nearDistanceTight, &
                           nearDistance, rand, nearBuf)
    else if (useHalfspace) then
      ! haveEval case (no patch data, but a real evaluate() exists)
      call sampleNearSurfEval(nTight, nLoose, nGauss, tgt, bboxMin, bboxMax, &
                               nearDistanceTight, nearDistance, rand, nearBuf)
    else
      batchSize = max(nNear * NEAR_BATCH_MULTIPLIER, MIN_BATCH_SIZE)
      nNearDone = 0
      do while (nNearDone < nNear)
        nBatch = min(batchSize, nNear - nNearDone) * NEAR_BATCH_MULTIPLIER
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
  end if

  ! ---- Write binary output ----
  if (nNear > 0) then
    call writeBinary(trim(outputFile), nVol, nNear, volBuf, nearBuf, bboxMin, bboxMax)
  else
    call writeBinary(trim(outputFile), nVol, 0, volBuf, volBuf(:, 1:0), bboxMin, bboxMax)
  end if

  ! ---- Clean up ----
  if (useUnion) then
    do i = 1, size(surfs)
      call surfs(i) % ptr % kill()
      deallocate(surfs(i) % ptr)
    end do
  else
    call surf % kill()
    deallocate(surf)
  end if
  deallocate(surfs)
  deallocate(patches)
  deallocate(volBuf)
  if (allocated(nearBuf)) deallocate(nearBuf)
  call input % kill()
  call surfDict % kill()
  call surfsDict % kill()

  print '(A,I0,A)', 'sdfSampler: wrote ', nVol + nNear, ' records to '//trim(outputFile)
  print '(A,I0,A,I0)', '  volumetric: ', nVol, '   near-surface: ', nNear
  if (nNear > 0) then
    if (useHalfspace) then
      print '(A,I0,A,I0,A,I0)', '  tight-linear: ', nTight, '   loose-linear: ', &
                                 nLoose, '   gaussian: ', nGauss
      print '(A,F8.4,A,F8.4,A,F8.4)', '  nearDistanceTight: ', nearDistanceTight, &
                                       '   nearDistance: ', nearDistance, &
                                       '   gaussianSigma: ', nearDistance * GAUSSIAN_SIGMA_FACTOR
    else
      print '(A,F8.4)', '  nearDelta: ', nearDelta
    end if
  end if

contains

  !!
  !! Halfspace test for a union of surfaces
  !!
  !! Matches unionCell's CSG union semantics: inside the union iff inside ANY
  !! constituent surface, i.e. outside iff outside EVERY constituent surface.
  !!
  !! Args:
  !!   surfs [in] -> Array of surfaces making up the union
  !!   r [in] -> Position
  !!   u [in] -> Normalised direction (norm2(u) = 1.0)
  !!
  !! Result:
  !!   .true. if r is outside every surface in the union, .false. if inside any
  !!
  function unionHalfspace(surfs, r, u) result(hs)
    type(surfPtr), dimension(:), intent(in) :: surfs
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)  :: hs
    integer(shortInt) :: i

    hs = .true.
    do i = 1, size(surfs)
      if (.not. surfs(i) % ptr % halfspace(r, u)) then
        hs = .false.
        return
      end if
    end do

  end function unionHalfspace

  !!
  !! Parse ctrlPts/weights/numPatches directly from a bezierVolume-format surface
  !! dictionary, mirroring bezierVolume_class's own init parsing (duplicated here
  !! since those arrays are private to the surface object, and sdfSampler only
  !! holds it via the polymorphic "surface" interface).
  !!
  !! Args:
  !!   dict [in] -> Surface sub-dictionary; must contain "numPatches" and "ctrlPts"
  !!   pd [out] -> Parsed patch data
  !!
  subroutine readPatchData(dict, pd)
    class(dictionary), intent(in) :: dict
    type(patchData), intent(out)  :: pd
    real(defReal), dimension(:), allocatable :: ctrlPtsList, weightsList
    integer(shortInt)                        :: m, p, j, k, l

    call dict % get(pd % numPatches, 'numPatches')
    call dict % get(ctrlPtsList,     'ctrlPts')

    allocate(pd % ctrlPts(pd % numPatches, 4, 4, 3))
    m = 1
    do p = 1, pd % numPatches
      do j = 1, 4
        do k = 1, 4
          do l = 1, 3
            pd % ctrlPts(p, j, k, l) = ctrlPtsList(m)
            m = m + 1
          end do
        end do
      end do
    end do

    allocate(pd % weights(pd % numPatches, 4, 4))
    if (dict % isPresent('weights')) then
      call dict % get(weightsList, 'weights')
      m = 1
      do p = 1, pd % numPatches
        do j = 1, 4
          do k = 1, 4
            pd % weights(p, j, k) = weightsList(m)
            m = m + 1
          end do
        end do
      end do
    else
      pd % weights = ONE
    end if

  end subroutine readPatchData

  !!
  !! Draw a standard-normal deviate via the Box-Muller transform
  !!
  !! Args:
  !!   rand [inout] -> Random number generator
  !!
  !! Result:
  !!   A single N(0,1) deviate
  !!
  function randNormal(rand) result(z)
    type(rng), intent(inout) :: rand
    real(defReal) :: z
    real(defReal) :: u1, u2

    u1 = rand % get()
    u2 = rand % get()
    u1 = max(u1, MIN_LOG_ARG)
    z = sqrt(-TWO * log(u1)) * cos(TWO_PI * u2)

  end function randNormal

  !!
  !! Draw a random point on a randomly chosen patch, with its unit normal
  !!
  !! The normal is estimated via finite differences in (u,v). Retries (up to
  !! maxAttempts) on a degenerate normal -- near-zero cross product, e.g.
  !! exactly at a fan-patch pole where all 4 corner control points coincide.
  !!
  !! Args:
  !!   patches [in] -> Raw patch data to sample from
  !!   maxAttempts [in] -> Safety cap on degenerate-normal retries
  !!   rand [inout] -> Random number generator
  !!   pC [out] -> Sampled surface point
  !!   nHat [out] -> Unit normal at pC
  !!
  !! Errors:
  !!   fatalError if maxAttempts is exceeded
  !!
  subroutine samplePatchPoint(patches, maxAttempts, rand, pC, nHat)
    type(patchData), dimension(:), intent(in) :: patches
    integer(shortInt), intent(in)             :: maxAttempts
    type(rng), intent(inout)                  :: rand
    real(defReal), dimension(3), intent(out)  :: pC, nHat
    integer(shortInt)           :: k, pIdx, attempt
    real(defReal)               :: uParam, vParam, uLo, uHi, vLo, vHi, normLen
    real(defReal), dimension(3) :: pU, pV
    character(100), parameter :: Here = 'samplePatchPoint (sdfSampler.f90)'

    attempt = 0
    do
      attempt = attempt + 1
      if (attempt > maxAttempts) then
        call fatalError(Here, 'Exceeded maxAttempts finding a non-degenerate surface &
                              &point (a wholly degenerate/zero-area patch would &
                              &trigger this).')
      end if

      k    = 1 + int(rand % get() * size(patches), shortInt)
      pIdx = 1 + int(rand % get() * patches(k) % numPatches, shortInt)
      uParam = rand % get()
      vParam = rand % get()

      uLo = max(uParam - FD_PARAM_EPS, ZERO)
      uHi = min(uParam + FD_PARAM_EPS, ONE)
      vLo = max(vParam - FD_PARAM_EPS, ZERO)
      vHi = min(vParam + FD_PARAM_EPS, ONE)

      pC = evalPatch(patches(k) % ctrlPts(pIdx,:,:,:), patches(k) % weights(pIdx,:,:), &
                     uParam, vParam)
      pU = evalPatch(patches(k) % ctrlPts(pIdx,:,:,:), patches(k) % weights(pIdx,:,:), &
                     uHi, vParam) &
         - evalPatch(patches(k) % ctrlPts(pIdx,:,:,:), patches(k) % weights(pIdx,:,:), &
                     uLo, vParam)
      pV = evalPatch(patches(k) % ctrlPts(pIdx,:,:,:), patches(k) % weights(pIdx,:,:), &
                     uParam, vHi) &
         - evalPatch(patches(k) % ctrlPts(pIdx,:,:,:), patches(k) % weights(pIdx,:,:), &
                     uParam, vLo)

      nHat = crossProduct(pU, pV)
      normLen = norm2(nHat)
      if (normLen >= MIN_NORMAL_LENGTH) exit
    end do
    nHat = nHat / normLen

  end subroutine samplePatchPoint

  !!
  !! Fill one near-surface band into nearBuf by sampling directly on the
  !! surface (random patch + random (u,v)) and offsetting along the normal,
  !! rather than by rejection sampling.
  !!
  !! Each accepted surface point yields two near-surface samples for the cost
  !! of two halfspace()/unionHalfspace() calls; the offset for each is drawn
  !! independently, Uniform(-dist,dist) if isGaussian is false or N(0,dist) if
  !! true. Labels always come from the FULL geometry (tgt), not just the
  !! sampled patch's own surface, so union overlap regions (e.g. the teapot
  !! handle re-entering the body) are still labelled correctly.
  !!
  !! Args:
  !!   nBand [in] -> Number of samples to draw for this band
  !!   isGaussian [in] -> True for a Gaussian(0,dist) offset, false for a
  !!     Uniform(-dist,dist) offset
  !!   dist [in] -> Gaussian standard deviation, or uniform half-width
  !!   tgt [in] -> Surface(s) to label the offset samples against
  !!   patches [in] -> Raw patch data to sample points from
  !!   rand [inout] -> Random number generator
  !!   nearBuf [inout] -> Buffer to fill with (x, y, z, +/-1 label)
  !!   startIdx [inout] -> Index of the last-filled column in nearBuf on entry;
  !!     advanced by nBand on exit
  !!
  subroutine fillBand(nBand, isGaussian, dist, tgt, patches, rand, nearBuf, startIdx)
    integer(shortInt), intent(in)                :: nBand
    logical(defBool), intent(in)                 :: isGaussian
    real(defReal), intent(in)                    :: dist
    type(sampleTarget), intent(in)               :: tgt
    type(patchData), dimension(:), intent(in)    :: patches
    type(rng), intent(inout)                     :: rand
    real(defReal), dimension(:,:), intent(inout) :: nearBuf
    integer(shortInt), intent(inout)             :: startIdx
    integer(shortInt)           :: nDone, maxAttempts
    real(defReal)               :: offset1, offset2
    real(defReal), dimension(3) :: pC, nHat, r1, r2, uLocal
    logical(defBool)            :: hs1, hs2

    if (nBand == 0) return

    uLocal = [ONE, ZERO, ZERO]
    maxAttempts = max(nBand * DEGENERACY_RETRY_MULTIPLIER, MIN_DEGENERACY_RETRIES)
    nDone = 0

    do while (nDone < nBand)
      call samplePatchPoint(patches, maxAttempts, rand, pC, nHat)

      if (isGaussian) then
        offset1 = randNormal(rand) * dist
        offset2 = randNormal(rand) * dist
      else
        offset1 = (TWO * rand % get() - ONE) * dist
        offset2 = (TWO * rand % get() - ONE) * dist
      end if

      r1 = pC + nHat * offset1
      r2 = pC + nHat * offset2

      if (tgt % useUnion) then
        hs1 = unionHalfspace(tgt % surfs, r1, uLocal)
        hs2 = unionHalfspace(tgt % surfs, r2, uLocal)
      else
        hs1 = tgt % surf % halfspace(r1, uLocal)
        hs2 = tgt % surf % halfspace(r2, uLocal)
      end if

      nDone = nDone + 1
      startIdx = startIdx + 1
      nearBuf(1:3, startIdx) = r1
      nearBuf(4,   startIdx) = merge(-ONE, ONE, hs1)
      if (nDone == nBand) exit

      nDone = nDone + 1
      startIdx = startIdx + 1
      nearBuf(1:3, startIdx) = r2
      nearBuf(4,   startIdx) = merge(-ONE, ONE, hs2)
    end do

  end subroutine fillBand

  !!
  !! Draw near-surface samples for halfspace-mode bezierVolume geometries
  !! across three bands -- tight linear, loose linear, Gaussian -- mirroring
  !! scripts/neuralSurface/sampler.py:generate_sphere_sdf's three-part
  !! near-surface strategy (same fractions and distance parameters, applied to
  !! Bezier patches via direct on-surface sampling instead of rejection
  !! sampling against an analytic sphere formula).
  !!
  !! Args:
  !!   nTight [in] -> Number of tight-linear-band samples
  !!   nLoose [in] -> Number of loose-linear-band samples
  !!   nGauss [in] -> Number of Gaussian-band samples
  !!   tgt [in] -> Surface(s) to label the samples against
  !!   patches [in] -> Raw patch data to sample points from
  !!   distTight [in] -> Tight-band uniform half-width
  !!   distLoose [in] -> Loose-band uniform half-width; the Gaussian band's
  !!     standard deviation is GAUSSIAN_SIGMA_FACTOR * distLoose
  !!   rand [inout] -> Random number generator
  !!   nearBuf [out] -> (4, nTight+nLoose+nGauss) buffer to fill with
  !!     (x, y, z, +/-1 label)
  !!
  subroutine sampleNearSurf(nTight, nLoose, nGauss, tgt, patches, distTight, distLoose, &
                             rand, nearBuf)
    integer(shortInt), intent(in)                                      :: nTight, nLoose, nGauss
    type(sampleTarget), intent(in)                                     :: tgt
    type(patchData), dimension(:), intent(in)                          :: patches
    real(defReal), intent(in)                                          :: distTight, distLoose
    type(rng), intent(inout)                                           :: rand
    real(defReal), dimension(4, nTight + nLoose + nGauss), intent(out) :: nearBuf
    integer(shortInt) :: startIdx

    startIdx = 0
    call fillBand(nTight, .false., distTight, tgt, patches, rand, nearBuf, startIdx)
    call fillBand(nLoose, .false., distLoose, tgt, patches, rand, nearBuf, startIdx)
    call fillBand(nGauss, .true., distLoose * GAUSSIAN_SIGMA_FACTOR, tgt, patches, &
                  rand, nearBuf, startIdx)

  end subroutine sampleNearSurf

  !!
  !! Draw a random point ON a surface's zero level set via Newton iteration on
  !! evaluate(), with the unit normal from the same finite-difference gradient.
  !! Generic across ANY surface with a genuine evaluate() (sphere, cylinder,
  !! plane, box, ...) -- no patch/control-point data needed, unlike
  !! samplePatchPoint. Convergence: for a quadratic residual F(r) (every
  !! current SCONE analytic surface), the Newton step p <- p - F(p)*grad(F)/
  !! |grad(F)|^2 is exactly Heron's method for the implied distance, so it
  !! converges quadratically fast (a handful of iterations) regardless of
  !! starting point, as long as the gradient is non-degenerate somewhere along
  !! the way.
  !!
  !! For a union, a random member surface is chosen each attempt (uniformly);
  !! the projection/normal come from that member alone (its own evaluate()),
  !! but the eventual LABEL (in fillBandEval) still comes from the full union,
  !! exactly as for the patch-based method.
  !!
  !! Args:
  !!   tgt [in] -> Surface(s) to project onto
  !!   bboxMin, bboxMax [in] -> Sampling bbox, used to draw the initial guess
  !!   maxAttempts [in] -> Safety cap on degenerate/non-converged retries
  !!   rand [inout] -> Random number generator
  !!   pC [out] -> Sampled surface point
  !!   nHat [out] -> Unit normal at pC
  !!
  !! Errors:
  !!   fatalError if maxAttempts is exceeded
  !!
  subroutine sampleEvalPoint(tgt, bboxMin, bboxMax, maxAttempts, rand, pC, nHat)
    type(sampleTarget), intent(in)           :: tgt
    real(defReal), dimension(3), intent(in)  :: bboxMin, bboxMax
    integer(shortInt), intent(in)            :: maxAttempts
    type(rng), intent(inout)                 :: rand
    real(defReal), dimension(3), intent(out) :: pC, nHat
    class(surface), pointer     :: activeSurf
    integer(shortInt)           :: attempt, iter, k
    real(defReal)               :: F0, normSq, normLen
    real(defReal), dimension(3) :: grad, p
    real(defReal), parameter     :: FD_EPS         = 1.0E-6_defReal
    real(defReal), parameter     :: NEWTON_TOL     = 1.0E-9_defReal
    real(defReal), parameter     :: CONVERGED_TOL  = 1.0E-4_defReal
    integer(shortInt), parameter :: MAX_NEWTON_ITERS = 30
    character(100), parameter    :: Here = 'sampleEvalPoint (sdfSampler.f90)'

    attempt = 0
    do
      attempt = attempt + 1
      if (attempt > maxAttempts) then
        call fatalError(Here, 'Exceeded maxAttempts finding a surface point via &
                              &evaluate()-based Newton projection.')
      end if

      if (tgt % useUnion) then
        k = 1 + int(rand % get() * size(tgt % surfs), shortInt)
        activeSurf => tgt % surfs(k) % ptr
      else
        activeSurf => tgt % surf
      end if

      p(1) = bboxMin(1) + rand % get() * (bboxMax(1) - bboxMin(1))
      p(2) = bboxMin(2) + rand % get() * (bboxMax(2) - bboxMin(2))
      p(3) = bboxMin(3) + rand % get() * (bboxMax(3) - bboxMin(3))

      do iter = 1, MAX_NEWTON_ITERS
        F0 = activeSurf % evaluate(p)
        grad = evalGradient(activeSurf, p, FD_EPS)
        normSq = dot_product(grad, grad)
        if (normSq < MIN_NORMAL_LENGTH) exit
        p = p - F0 * grad / normSq
        if (abs(F0) < NEWTON_TOL) exit
      end do

      F0     = activeSurf % evaluate(p)
      grad   = evalGradient(activeSurf, p, FD_EPS)
      normLen = norm2(grad)
      if (normLen >= MIN_NORMAL_LENGTH .and. abs(F0) < CONVERGED_TOL) then
        pC   = p
        nHat = grad / normLen
        exit
      end if
    end do

  end subroutine sampleEvalPoint

  !!
  !! Central-difference gradient of a surface's evaluate() at a point
  !!
  !! Args:
  !!   activeSurf [in] -> Surface whose evaluate() is differentiated
  !!   p [in] -> Point at which to estimate the gradient
  !!   eps [in] -> Finite-difference step size
  !!
  !! Result:
  !!   Estimated gradient of evaluate() at p
  !!
  function evalGradient(activeSurf, p, eps) result(grad)
    class(surface), pointer, intent(in)     :: activeSurf
    real(defReal), dimension(3), intent(in) :: p
    real(defReal), intent(in)               :: eps
    real(defReal), dimension(3) :: grad
    integer(shortInt)           :: i
    real(defReal), dimension(3) :: pPlus, pMinus

    do i = 1, 3
      pPlus    = p
      pMinus   = p
      pPlus(i)  = pPlus(i)  + eps
      pMinus(i) = pMinus(i) - eps
      grad(i) = (activeSurf % evaluate(pPlus) - activeSurf % evaluate(pMinus)) / (TWO * eps)
    end do

  end function evalGradient

  !!
  !! Generic (evaluate()-based) analogue of fillBand for surfaces without
  !! patch data -- see sampleEvalPoint.
  !!
  !! Args:
  !!   nBand [in] -> Number of samples to draw for this band
  !!   isGaussian [in] -> True for a Gaussian(0,dist) offset, false for a
  !!     Uniform(-dist,dist) offset
  !!   dist [in] -> Gaussian standard deviation, or uniform half-width
  !!   tgt [in] -> Surface(s) to label the offset samples against
  !!   bboxMin, bboxMax [in] -> Sampling bbox, used to draw the Newton-projection
  !!     initial guess
  !!   rand [inout] -> Random number generator
  !!   nearBuf [inout] -> Buffer to fill with (x, y, z, +/-1 label)
  !!   startIdx [inout] -> Index of the last-filled column in nearBuf on entry;
  !!     advanced by nBand on exit
  !!
  subroutine fillBandEval(nBand, isGaussian, dist, tgt, bboxMin, bboxMax, rand, nearBuf, &
                           startIdx)
    integer(shortInt), intent(in)                :: nBand
    logical(defBool), intent(in)                 :: isGaussian
    real(defReal), intent(in)                    :: dist
    type(sampleTarget), intent(in)               :: tgt
    real(defReal), dimension(3), intent(in)      :: bboxMin, bboxMax
    type(rng), intent(inout)                     :: rand
    real(defReal), dimension(:,:), intent(inout) :: nearBuf
    integer(shortInt), intent(inout)             :: startIdx
    integer(shortInt)           :: nDone, maxAttempts
    real(defReal)               :: offset1, offset2
    real(defReal), dimension(3) :: pC, nHat, r1, r2, uLocal
    logical(defBool)            :: hs1, hs2

    if (nBand == 0) return

    uLocal = [ONE, ZERO, ZERO]
    maxAttempts = max(nBand * DEGENERACY_RETRY_MULTIPLIER, MIN_DEGENERACY_RETRIES)
    nDone = 0

    do while (nDone < nBand)
      call sampleEvalPoint(tgt, bboxMin, bboxMax, maxAttempts, rand, pC, nHat)

      if (isGaussian) then
        offset1 = randNormal(rand) * dist
        offset2 = randNormal(rand) * dist
      else
        offset1 = (TWO * rand % get() - ONE) * dist
        offset2 = (TWO * rand % get() - ONE) * dist
      end if

      r1 = pC + nHat * offset1
      r2 = pC + nHat * offset2

      if (tgt % useUnion) then
        hs1 = unionHalfspace(tgt % surfs, r1, uLocal)
        hs2 = unionHalfspace(tgt % surfs, r2, uLocal)
      else
        hs1 = tgt % surf % halfspace(r1, uLocal)
        hs2 = tgt % surf % halfspace(r2, uLocal)
      end if

      nDone = nDone + 1
      startIdx = startIdx + 1
      nearBuf(1:3, startIdx) = r1
      nearBuf(4,   startIdx) = merge(-ONE, ONE, hs1)
      if (nDone == nBand) exit

      nDone = nDone + 1
      startIdx = startIdx + 1
      nearBuf(1:3, startIdx) = r2
      nearBuf(4,   startIdx) = merge(-ONE, ONE, hs2)
    end do

  end subroutine fillBandEval

  !!
  !! Generic (evaluate()-based) analogue of sampleNearSurf for surfaces
  !! without patch data -- see sampleEvalPoint.
  !!
  !! Args:
  !!   nTight [in] -> Number of tight-linear-band samples
  !!   nLoose [in] -> Number of loose-linear-band samples
  !!   nGauss [in] -> Number of Gaussian-band samples
  !!   tgt [in] -> Surface(s) to label the samples against
  !!   bboxMin, bboxMax [in] -> Sampling bbox, used to draw the Newton-projection
  !!     initial guess
  !!   distTight [in] -> Tight-band uniform half-width
  !!   distLoose [in] -> Loose-band uniform half-width; the Gaussian band's
  !!     standard deviation is GAUSSIAN_SIGMA_FACTOR * distLoose
  !!   rand [inout] -> Random number generator
  !!   nearBuf [out] -> (4, nTight+nLoose+nGauss) buffer to fill with
  !!     (x, y, z, +/-1 label)
  !!
  subroutine sampleNearSurfEval(nTight, nLoose, nGauss, tgt, bboxMin, bboxMax, distTight, &
                                 distLoose, rand, nearBuf)
    integer(shortInt), intent(in)                                      :: nTight, nLoose, nGauss
    type(sampleTarget), intent(in)                                     :: tgt
    real(defReal), dimension(3), intent(in)                            :: bboxMin, bboxMax
    real(defReal), intent(in)                                          :: distTight, distLoose
    type(rng), intent(inout)                                           :: rand
    real(defReal), dimension(4, nTight + nLoose + nGauss), intent(out) :: nearBuf
    integer(shortInt) :: startIdx

    startIdx = 0
    call fillBandEval(nTight, .false., distTight, tgt, bboxMin, bboxMax, rand, nearBuf, startIdx)
    call fillBandEval(nLoose, .false., distLoose, tgt, bboxMin, bboxMax, rand, nearBuf, startIdx)
    call fillBandEval(nGauss, .true., distLoose * GAUSSIAN_SIGMA_FACTOR, tgt, bboxMin, bboxMax, &
                       rand, nearBuf, startIdx)

  end subroutine sampleNearSurfEval

  !!
  !! Write the sampled points to a binary file readable by
  !! scripts/neuralSurface/sampler.py:load_scone_binary()
  !!
  !! Args:
  !!   filename [in] -> Path of the output file
  !!   nVol [in] -> Number of volumetric-sample records
  !!   nNear [in] -> Number of near-surface-sample records
  !!   volBuf [in] -> (4, nVol) volumetric samples: x, y, z, sdf/label
  !!   nearBuf [in] -> (4, nNear) near-surface samples: x, y, z, sdf/label
  !!   bboxMin, bboxMax [in] -> Sampling bbox, written into the file header
  !!
  !! Errors:
  !!   fatalError if the output file cannot be opened
  !!
  subroutine writeBinary(filename, nVol, nNear, volBuf, nearBuf, bboxMin, bboxMax)
    character(*), intent(in)                :: filename
    integer(shortInt), intent(in)           :: nVol, nNear
    real(defReal), intent(in)               :: volBuf(4, nVol)
    real(defReal), intent(in)               :: nearBuf(4, nNear)
    real(defReal), dimension(3), intent(in) :: bboxMin, bboxMax
    integer(shortInt) :: i, unit
    integer(int32)    :: count32
    real(real64)      :: bbox8(6), rec8(4)
    character(100), parameter :: Here = 'writeBinary (sdfSampler.f90)'

    count32 = int(nVol + nNear, int32)
    bbox8(1:3) = real(bboxMin, real64)
    bbox8(4:6) = real(bboxMax, real64)

    open(newunit=unit, file=filename, status='replace', &
         access='stream', form='unformatted', iostat=i)
    if (i /= 0) call fatalError(Here, 'Cannot open output file: '//trim(filename))

    write(unit) count32
    write(unit) bbox8
    do i = 1, nVol
      rec8 = real(volBuf(:, i), real64)
      write(unit) rec8
    end do
    do i = 1, nNear
      rec8 = real(nearBuf(:, i), real64)
      write(unit) rec8
    end do
    close(unit)

  end subroutine writeBinary

end program sdfSampler
