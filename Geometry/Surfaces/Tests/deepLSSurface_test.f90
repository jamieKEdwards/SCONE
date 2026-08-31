module deepLSSurface_test

  use numPrecision
  use universalVariables, only : INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use trainedMLP_class,   only : trainedMLP, ACTIVATION_LEAKYRELU
  use deepLSSurface_class, only : deepLSSurface
  use funit

  implicit none

  !!
  !! On-disk v3 DeepLS weight fixture (see scripts/neuralSurface/export_deepls.py
  !! for the byte layout; the reader is deepLSSurface_class % loadWeights).
  !!
  !! Grid: nvox = [1,1,3], gridOrigin = [-1.5,-1.5,-1.5], voxelSize = [3,3,1]
  !!   -> world domain x,y in [-1.5,1.5], z split into three unit-thick bands:
  !!        voxel (1,1,1): z in [-1.5,-0.5]   status CONST_INSIDE  -> F = -1
  !!        voxel (1,1,2): z in [-0.5, 0.5]   status HAS_MLP
  !!        voxel (1,1,3): z in [ 0.5, 1.5]   status CONST_OUTSIDE -> F = +1
  !!
  !! Shared decoder (inputDim = latentDim + 3 = 4, hidden 2, 2 layers, LeakyReLU):
  !!   Layer 1: W(1,:) = [0,1,0,0], W(2,:) = [0,0,1,0]  (selects local xn, yn;
  !!            zeroes the latent column and the local zn column)
  !!   Layer 2: W(1,:) = [1,1]
  !!   => decoder([lat, xn, yn, zn]) = tanh( LeakyReLU(xn) + LeakyReLU(yn) )
  !!
  !! Active voxel (1,1,2): local bbox [-1,1]^3, latent = [0.2].
  !! For world point [0.5,0.5,0.0] (rs = the same at geomScale 1): xn,yn = 0.5,
  !!   so F = tanh(1.0).
  !!
  character(*), parameter :: FIXTURE    = '/tmp/scone_deepLSSurface_test.bin'
  real(defReal), parameter :: GEOM_SCALE = 2.0_defReal
  real(defReal), parameter :: TOL        = 1.0e-12_defReal

  type(deepLSSurface) :: surf        ! geomScale = 1
  type(deepLSSurface) :: surfScaled  ! geomScale = GEOM_SCALE

contains

  !!
  !! Hand-write the v3 binary fixture (deepLSWeightIO has no writer)
  !!
  subroutine writeFixture(path)
    character(*), intent(in)      :: path
    integer                       :: u, l
    integer(shortInt), parameter  :: MAGIC   = int(z'4E534446', shortInt)
    integer(shortInt), parameter  :: VERSION = 3_shortInt
    integer(shortInt), parameter  :: LAT = 1_shortInt, HID = 2_shortInt, NL = 2_shortInt
    integer(shortInt)             :: statusMap(1, 1, 3)
    integer(shortInt)             :: inD, outD
    type(trainedMLP)              :: dec
    real(defReal)                 :: testOut

    call dec % init(LAT + 3_shortInt, HID, NL, ACTIVATION_LEAKYRELU, 0.01_defReal, ONE, &
                    [ZERO, ZERO, ZERO], [ONE, ONE, ONE])
    dec % weights(:, :, :) = ZERO
    dec % weights(1, 2, 1) = ONE     ! hidden 1 <- local xn
    dec % weights(2, 3, 1) = ONE     ! hidden 2 <- local yn
    dec % weights(1, 1, 2) = ONE     ! output <- hidden 1
    dec % weights(1, 2, 2) = ONE     ! output <- hidden 2

    statusMap(1, 1, 1) = 0_shortInt  ! CONST_INSIDE
    statusMap(1, 1, 2) = 2_shortInt  ! HAS_MLP
    statusMap(1, 1, 3) = 1_shortInt  ! CONST_OUTSIDE

    testOut = tanh(ONE)

    open(newunit=u, file=path, access='stream', form='unformatted', &
         status='replace', action='write')

    ! Header
    write(u) MAGIC
    write(u) VERSION
    write(u) LAT
    write(u) HID
    write(u) NL
    write(u) ACTIVATION_LEAKYRELU
    write(u) 0.01_defReal                                       ! leaky_alpha
    write(u) 1_shortInt, 1_shortInt, 3_shortInt                 ! nvox
    write(u) [-1.5_defReal, -1.5_defReal, -1.5_defReal]         ! grid_origin
    write(u) [3.0_defReal, 3.0_defReal, 1.0_defReal]            ! voxel_size
    write(u) 1_shortInt                                         ! n_active

    ! Voxel map (Fortran column-major, ix fastest)
    write(u) statusMap

    ! Shared decoder, layer by layer (weight (out,in) col-major, then bias(out))
    do l = 1, int(NL)
      call dec % layerDims(int(l, shortInt), inD, outD)
      write(u) dec % weights(1:outD, 1:inD, l)
      write(u) dec % biases(1:outD, l)
    end do

    ! Active voxel latent block
    write(u) 1_shortInt, 1_shortInt, 2_shortInt                 ! voxel_index (1-based)
    write(u) [-ONE, -ONE, -ONE]                                 ! bbox_min
    write(u) [ONE, ONE, ONE]                                    ! bbox_max
    write(u) [0.2_defReal]                                      ! latent (LAT = 1)

    ! Validation test vector (first active voxel)
    write(u) [0.5_defReal, 0.5_defReal, ZERO]                   ! test_input
    write(u) testOut                                            ! test_output

    close(u)
    call dec % kill()

  end subroutine writeFixture

@Before
  subroutine setUp()
    type(dictionary) :: dict

    call writeFixture(FIXTURE)

    call charToDict(dict, "id 1; weightFile " // FIXTURE // ";")
    call surf % init(dict)
    call dict % kill()

    call charToDict(dict, "id 2; weightFile " // FIXTURE // "; geometricScale 2.0;")
    call surfScaled % init(dict)
    call dict % kill()

  end subroutine setUp

@After
  subroutine tearDown()
    integer(shortInt) :: u, stat
    logical(defBool)  :: exists

    call surf % kill()
    call surfScaled % kill()

    ! Remove the fixture so a stale file can never mask a reader/writer regression
    inquire(file=FIXTURE, exist=exists)
    if (.not. exists) return
    open(newunit=u, file=FIXTURE, status='old', iostat=stat)
    if (stat == 0) close(u, status='delete')

  end subroutine tearDown

  !!
  !! Surface type name
  !!
@Test
  subroutine testMyType()

    @assertEqual('deepLSSurface', surf % myType())

  end subroutine testMyType

  !!
  !! Bounding box (default geomScale) is the voxel-grid extent
  !!
@Test
  subroutine testBoundingBoxUnitScale()
    real(defReal), dimension(6) :: aabb

    aabb = surf % boundingBox()
    @assertEqual([-1.5_defReal, -1.5_defReal, -1.5_defReal, &
                   1.5_defReal,  1.5_defReal,  1.5_defReal], aabb, TOL)

  end subroutine testBoundingBoxUnitScale

  !!
  !! Bounding box scales with geomScale
  !!
@Test
  subroutine testBoundingBoxWithGeomScale()
    real(defReal), dimension(6) :: aabb

    aabb = surfScaled % boundingBox()
    @assertEqual([-3.0_defReal, -3.0_defReal, -3.0_defReal, &
                   3.0_defReal,  3.0_defReal,  3.0_defReal], aabb, TOL)

  end subroutine testBoundingBoxWithGeomScale

  !!
  !! evaluate() branch: CONST_INSIDE voxel returns -1
  !!
@Test
  subroutine testEvaluateConstInsideVoxel()

    ! z = -1.0 lands in voxel (1,1,1)
    @assertEqual(-ONE, surf % evaluate([ZERO, ZERO, -ONE]), TOL)

  end subroutine testEvaluateConstInsideVoxel

  !!
  !! evaluate() branch: HAS_MLP voxel runs the shared decoder
  !!
@Test
  subroutine testEvaluateHasMlpVoxel()

    ! [0.5,0.5,0] lands in voxel (1,1,2); decoder -> tanh(0.5 + 0.5)
    @assertEqual(tanh(ONE), surf % evaluate([0.5_defReal, 0.5_defReal, ZERO]), TOL)

  end subroutine testEvaluateHasMlpVoxel

  !!
  !! evaluate() branch: CONST_OUTSIDE voxel returns +1
  !!
@Test
  subroutine testEvaluateConstOutsideVoxel()

    ! z = 1.0 lands in voxel (1,1,3)
    @assertEqual(ONE, surf % evaluate([ZERO, ZERO, ONE]), TOL)

  end subroutine testEvaluateConstOutsideVoxel

  !!
  !! evaluate() branch: outside the voxel grid returns +1
  !!
@Test
  subroutine testEvaluateOutsideGrid()

    ! z = 5.0 is above the grid (z_max = 1.5)
    @assertEqual(ONE, surf % evaluate([ZERO, ZERO, 5.0_defReal]), TOL)

  end subroutine testEvaluateOutsideGrid

  !!
  !! evaluate() with geomScale: coords divided by s before the grid lookup,
  !! result multiplied back by s
  !!
@Test
  subroutine testEvaluateWithGeomScale()

    ! r = [1,1,0] -> rs = [0.5,0.5,0] -> HAS_MLP voxel -> tanh(1), * s
    @assertEqual(GEOM_SCALE * tanh(ONE), &
                 surfScaled % evaluate([ONE, ONE, ZERO]), TOL)

  end subroutine testEvaluateWithGeomScale

  !!
  !! distance() is a stub returning INF (delta-tracking-only surface)
  !!
@Test
  subroutine testDistanceReturnsInf()

    @assertGreaterThanOrEqual(surf % distance([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO]), INF)

  end subroutine testDistanceReturnsInf

  !!
  !! going(): forward finite-difference step resolves the halfspace
  !!
@Test
  subroutine testGoingPositive()

    ! From [0,0,0] (F = 0 in the HAS_MLP voxel) stepping +x lands in F > 0
    @assertTrue(surf % going([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO]))

  end subroutine testGoingPositive

@Test
  subroutine testGoingNegative()

    ! Stepping -x from [0,0,0] lands in F < 0 (LeakyReLU negative slope)
    @assertFalse(surf % going([ZERO, ZERO, ZERO], [-ONE, ZERO, ZERO]))

  end subroutine testGoingNegative

  !!
  !! halfspace(): sign of evaluate() through the const voxels
  !!
@Test
  subroutine testHalfspacePositive()

    @assertTrue(surf % halfspace([ZERO, ZERO, ONE], [ONE, ZERO, ZERO]))

  end subroutine testHalfspacePositive

@Test
  subroutine testHalfspaceNegative()

    @assertFalse(surf % halfspace([ZERO, ZERO, -ONE], [ONE, ZERO, ZERO]))

  end subroutine testHalfspaceNegative

end module deepLSSurface_test
