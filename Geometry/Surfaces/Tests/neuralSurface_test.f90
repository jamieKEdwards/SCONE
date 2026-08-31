module neuralSurface_test

  use numPrecision
  use universalVariables, only : INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use trainedMLP_class,   only : trainedMLP, ACTIVATION_LEAKYRELU
  use neuralSurface_class, only : neuralSurface
  use funit

  implicit none

  !!
  !! Reference MLP used to build the on-disk weight fixtures:
  !!   Architecture: inputDim=3, hiddenDim=2, numLayers=2, LeakyReLU(0.01)
  !!   Layer 1: W = [[1,0,0],[0,1,0]], b = [0,0]   (3 -> 2)
  !!   Layer 2: W = [[1,1]],           b = [0]     (2 -> 1)
  !!   bbox: [-1,-1,-1] .. [1,1,1], sdfScale = 1
  !!
  !! For a point p (inside the bbox, so x_norm = p):
  !!   F(p) = tanh( LeakyReLU(p_x) + LeakyReLU(p_y) )
  !! e.g. F([0.5,0.5,0.5]) = tanh(1.0), F([0,0,0]) = 0.
  !!
  character(*), parameter :: FIX_UNIT   = '/tmp/scone_neuralSurface_test_unit.bin'
  character(*), parameter :: FIX_SCALED = '/tmp/scone_neuralSurface_test_scaled.bin'

  ! geomScale used by the scaled instance
  real(defReal), parameter :: GEOM_SCALE = 2.0_defReal
  real(defReal), parameter :: TOL        = 1.0e-12_defReal

  type(neuralSurface) :: surf        ! geomScale = 1 (dict omits geometricScale)
  type(neuralSurface) :: surfScaled  ! geomScale = GEOM_SCALE

contains

  !!
  !! Write the reference MLP to a binary weight fixture on disk
  !!
  subroutine writeFixture(path)
    character(*), intent(in)    :: path
    type(trainedMLP)            :: mlp
    real(defReal), dimension(3) :: testInput
    real(defReal)               :: testOutput

    call mlp % init(inputDim       = 3,                    &
                    hiddenDim      = 2,                    &
                    numLayers      = 2,                    &
                    activationType = ACTIVATION_LEAKYRELU, &
                    leakyAlpha     = 0.01_defReal,         &
                    sdfScale       = ONE,                  &
                    bboxMin        = [-ONE, -ONE, -ONE],   &
                    bboxMax        = [ ONE,  ONE,  ONE])

    mlp % weights(1, 1, 1) = ONE;  mlp % weights(1, 2, 1) = ZERO; mlp % weights(1, 3, 1) = ZERO
    mlp % weights(2, 1, 1) = ZERO; mlp % weights(2, 2, 1) = ONE;  mlp % weights(2, 3, 1) = ZERO
    mlp % weights(1, 1, 2) = ONE;  mlp % weights(1, 2, 2) = ONE

    testInput  = [0.5_defReal, 0.5_defReal, 0.5_defReal]
    testOutput = tanh(ONE)
    call mlp % dump(path, testInput, testOutput)
    call mlp % kill()

  end subroutine writeFixture

  !!
  !! Build both surface instances from freshly written fixtures
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    call writeFixture(FIX_UNIT)
    call writeFixture(FIX_SCALED)

    call charToDict(dict, "id 1; weightFile " // FIX_UNIT // ";")
    call surf % init(dict)
    call dict % kill()

    call charToDict(dict, "id 2; weightFile " // FIX_SCALED // "; geometricScale 2.0;")
    call surfScaled % init(dict)
    call dict % kill()

  end subroutine setUp

@After
  subroutine tearDown()

    call surf % kill()
    call surfScaled % kill()

    ! Remove the fixtures so a stale file can never mask a writer regression
    call deleteFile(FIX_UNIT)
    call deleteFile(FIX_SCALED)

  end subroutine tearDown

  subroutine deleteFile(path)
    character(*), intent(in)      :: path
    integer(shortInt)             :: unit, stat
    logical(defBool)              :: exists

    inquire(file=path, exist=exists)
    if (.not. exists) return
    open(newunit=unit, file=path, status='old', iostat=stat)
    if (stat == 0) close(unit, status='delete')

  end subroutine deleteFile

  !!
  !! Surface type name
  !!
@Test
  subroutine testMyType()

    @assertEqual('neuralSurface', surf % myType())

  end subroutine testMyType

  !!
  !! Bounding box with the default geomScale = 1 matches the MLP's own bbox
  !!
@Test
  subroutine testBoundingBoxUnitScale()
    real(defReal), dimension(6) :: aabb

    aabb = surf % boundingBox()
    @assertEqual([-ONE, -ONE, -ONE, ONE, ONE, ONE], aabb, TOL)

  end subroutine testBoundingBoxUnitScale

  !!
  !! Bounding box scales with geomScale
  !!
@Test
  subroutine testBoundingBoxWithGeomScale()
    real(defReal), dimension(6) :: aabb

    aabb = surfScaled % boundingBox()
    @assertEqual([-TWO, -TWO, -TWO, TWO, TWO, TWO], aabb, TOL)

  end subroutine testBoundingBoxWithGeomScale

  !!
  !! evaluate() delegates to the MLP forward pass (geomScale = 1)
  !!
@Test
  subroutine testEvaluateUnitScale()

    @assertEqual(tanh(ONE), surf % evaluate([0.5_defReal, 0.5_defReal, 0.5_defReal]), TOL)

  end subroutine testEvaluateUnitScale

  !!
  !! evaluate() at the bbox centre is exactly zero
  !!
@Test
  subroutine testEvaluateAtBboxCentre()

    @assertEqual(ZERO, surf % evaluate([ZERO, ZERO, ZERO]), TOL)

  end subroutine testEvaluateAtBboxCentre

  !!
  !! evaluate() with geomScale: F(r) = mlp % evaluate(r / s) * s
  !!
  !! For r = [1,1,1] and s = 2 the MLP sees [0.5,0.5,0.5] (-> tanh(1)),
  !! scaled back up by 2.
  !!
@Test
  subroutine testEvaluateWithGeomScale()

    @assertEqual(TWO * tanh(ONE), &
                 surfScaled % evaluate([ONE, ONE, ONE]), TOL)

  end subroutine testEvaluateWithGeomScale

  !!
  !! distance() is a stub returning INF (delta-tracking-only surface)
  !!
@Test
  subroutine testDistanceReturnsInf()

    @assertGreaterThanOrEqual(surf % distance([0.5_defReal, ZERO, ZERO], [ONE, ZERO, ZERO]), INF)

  end subroutine testDistanceReturnsInf

  !!
  !! going(): forward finite-difference step along u resolves the halfspace
  !!
@Test
  subroutine testGoingPositive()

    ! From the centre (F = 0) stepping along +x lands in F > 0
    @assertTrue(surf % going([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO]))

  end subroutine testGoingPositive

@Test
  subroutine testGoingNegative()

    ! Stepping along -x from the centre lands in F < 0 (LeakyReLU negative slope)
    @assertFalse(surf % going([ZERO, ZERO, ZERO], [-ONE, ZERO, ZERO]))

  end subroutine testGoingNegative

  !!
  !! halfspace(): sign of evaluate() away from the surface tolerance band
  !!
@Test
  subroutine testHalfspacePositive()

    @assertTrue(surf % halfspace([0.5_defReal, 0.5_defReal, 0.5_defReal], [ONE, ZERO, ZERO]))

  end subroutine testHalfspacePositive

@Test
  subroutine testHalfspaceNegative()

    @assertFalse(surf % halfspace([-0.5_defReal, -0.5_defReal, -0.5_defReal], [ONE, ZERO, ZERO]))

  end subroutine testHalfspaceNegative

end module neuralSurface_test
