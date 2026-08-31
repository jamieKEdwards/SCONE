module trainedMLP_test

  use numPrecision
  use trainedMLP_class, only : trainedMLP, ACTIVATION_LEAKYRELU, ACTIVATION_RELU, &
                                ACTIVATION_TANH, MLP_MAX_DIM
  use funit

  implicit none

  !!
  !! Minimal 2-layer MLP for tests:
  !!   Architecture: inputDim=3, hiddenDim=2, numLayers=2
  !!   Layer 1: W = [[1,0,0],[0,1,0]], b = [0,0]  (3->2)
  !!   Layer 2: W = [[1,1]],           b = [0]    (2->1)
  !!   Activation: LeakyReLU (alpha=0.01), sdfScale=1.0
  !!   bbox: [-1,-1,-1] to [1,1,1]
  !!
  !! For input [0.5, 0.5, 0.5]:
  !!   normalised: [0.5, 0.5, 0.5]  (bbox is [-1,1] so x_norm = x)
  !!   h1 = W1 * [0.5,0.5,0.5] + b1 = [0.5, 0.5]
  !!   after LeakyReLU: [0.5, 0.5]
  !!   h2 = W2 * [0.5, 0.5] + b2 = [1.0]
  !!   sdf = tanh(1.0) * 1.0 = 0.76159...
  !!
  type(trainedMLP), save :: mlp

contains

  !!
  !! Build the test MLP before each test
  !!
@Before
  subroutine setUp()
    real(defReal), dimension(3) :: bboxMin, bboxMax

    bboxMin = [-ONE, -ONE, -ONE]
    bboxMax = [ ONE,  ONE,  ONE]

    call mlp % init(inputDim       = 3,                   &
                    hiddenDim      = 2,                   &
                    numLayers      = 2,                   &
                    activationType = ACTIVATION_LEAKYRELU,&
                    leakyAlpha     = 0.01_defReal,        &
                    sdfScale       = ONE,                 &
                    bboxMin        = bboxMin,             &
                    bboxMax        = bboxMax)

    ! Layer 1 weights: W(1:2, 1:3, 1)
    ! Row 1: [1, 0, 0]
    ! Row 2: [0, 1, 0]
    mlp % weights(1, 1, 1) = ONE;  mlp % weights(1, 2, 1) = ZERO; mlp % weights(1, 3, 1) = ZERO
    mlp % weights(2, 1, 1) = ZERO; mlp % weights(2, 2, 1) = ONE;  mlp % weights(2, 3, 1) = ZERO
    mlp % biases(1:2, 1) = ZERO

    ! Layer 2 weights: W(1:1, 1:2, 2)
    ! Row 1: [1, 1]
    mlp % weights(1, 1, 2) = ONE; mlp % weights(1, 2, 2) = ONE
    mlp % biases(1, 2) = ZERO

  end subroutine setUp

  !!
  !! Clean up after each test
  !!
@After
  subroutine cleanUp()
    call mlp % kill()
  end subroutine cleanUp

  !!
  !! Test that init() allocates correctly and sets isInit flag
  !!
@Test
  subroutine testInit()

    @assertTrue(mlp % isInit)
    @assertEqual(3, mlp % inputDim)
    @assertEqual(2, mlp % hiddenDim)
    @assertEqual(2, mlp % numLayers)
    @assertTrue(allocated(mlp % weights))
    @assertTrue(allocated(mlp % biases))

  end subroutine testInit

  !!
  !! Test forward pass with positive inputs (no LeakyReLU activation needed)
  !!
@Test
  subroutine testForwardPositive()
    real(defReal), dimension(3) :: point
    real(defReal)               :: sdf, expected
    real(defReal), parameter    :: TOL = 1.0e-12_defReal

    point = [0.5_defReal, 0.5_defReal, 0.5_defReal]
    sdf   = mlp % evaluate(point)

    ! Layer 1: [0.5, 0.5], after LeakyReLU: [0.5, 0.5]
    ! Layer 2: [1.0], sdf = tanh(1.0)
    expected = tanh(ONE)
    @assertEqual(expected, sdf, TOL)

  end subroutine testForwardPositive

  !!
  !! Test forward pass with negative inputs — exercises LeakyReLU negative slope
  !!
@Test
  subroutine testForwardNegativeLeakyReLU()
    real(defReal), dimension(3) :: point
    real(defReal)               :: sdf, expected
    real(defReal), parameter    :: TOL = 1.0e-12_defReal

    ! Input outside bbox: normalised x_norm = 2*(x - (-1))/(2) - 1 = x
    ! For x = -0.5: x_norm = -0.5
    point = [-0.5_defReal, -0.5_defReal, -0.5_defReal]
    sdf   = mlp % evaluate(point)

    ! Layer 1: W*[-0.5,-0.5,-0.5] = [-0.5, -0.5]
    ! After LeakyReLU(alpha=0.01): [-0.005, -0.005]
    ! Layer 2: W*[-0.005,-0.005] = [-0.01]
    ! sdf = tanh(-0.01) * 1.0
    expected = tanh(-0.01_defReal)
    @assertEqual(expected, sdf, TOL)

  end subroutine testForwardNegativeLeakyReLU

  !!
  !! Test coordinate normalisation: centre of bbox maps to [0,0,0]
  !!
@Test
  subroutine testNormalisation()
    real(defReal), dimension(3) :: point
    real(defReal)               :: sdf, expected
    real(defReal), parameter    :: TOL = 1.0e-12_defReal

    ! bbox = [-1,1]^3, centre = [0,0,0]
    ! x_norm = 2*(0-(-1))/(1-(-1)) - 1 = 2*1/2 - 1 = 0
    ! Layer 1: W*[0,0,0] = [0,0], after LeakyReLU = [0,0]
    ! Layer 2: [0], sdf = tanh(0) = 0
    point    = [ZERO, ZERO, ZERO]
    sdf      = mlp % evaluate(point)
    expected = ZERO
    @assertEqual(expected, sdf, TOL)

  end subroutine testNormalisation

  !!
  !! Test tanh output is bounded within [-sdfScale, sdfScale]
  !!
@Test
  subroutine testOutputBounded()
    real(defReal), dimension(3) :: point
    real(defReal)               :: sdf

    ! Use a large input to saturate the network
    point = [100.0_defReal, 100.0_defReal, 100.0_defReal]
    sdf   = mlp % evaluate(point)

    ! sdfScale = 1.0, so output must be in [-1, 1] (tanh saturates to exactly +-1.0)
    @assertTrue(sdf >= -ONE)
    @assertTrue(sdf <=  ONE)

    point = [-100.0_defReal, -100.0_defReal, -100.0_defReal]
    sdf   = mlp % evaluate(point)
    @assertTrue(sdf >= -ONE)
    @assertTrue(sdf <=  ONE)

  end subroutine testOutputBounded

  !!
  !! Test kill() deallocates and resets state
  !!
@Test
  subroutine testKill()

    call mlp % kill()
    @assertFalse(mlp % isInit)
    @assertFalse(allocated(mlp % weights))
    @assertFalse(allocated(mlp % biases))
    @assertEqual(0, mlp % inputDim)

    ! Re-init for cleanUp
    call setUp()

  end subroutine testKill

  !!
  !! Test with sdfScale != 1.0 — output should scale accordingly
  !!
@Test
  subroutine testSdfScale()
    type(trainedMLP)            :: mlp2
    real(defReal), dimension(3) :: bboxMin, bboxMax, point
    real(defReal)               :: sdf_scaled, sdf_unit
    real(defReal), parameter    :: SCALE = 5.0_defReal
    real(defReal), parameter    :: TOL   = 1.0e-12_defReal

    bboxMin = [-ONE, -ONE, -ONE]
    bboxMax = [ ONE,  ONE,  ONE]
    call mlp2 % init(3, 2, 2, ACTIVATION_LEAKYRELU, 0.01_defReal, SCALE, bboxMin, bboxMax)

    ! Same weights as module mlp
    mlp2 % weights = mlp % weights
    mlp2 % biases  = mlp % biases

    point       = [0.5_defReal, 0.5_defReal, 0.5_defReal]
    sdf_unit    = mlp  % evaluate(point)
    sdf_scaled  = mlp2 % evaluate(point)

    ! sdf_scaled should equal SCALE * sdf_unit (since tanh argument is the same)
    @assertEqual(SCALE * sdf_unit, sdf_scaled, TOL)

    call mlp2 % kill()

  end subroutine testSdfScale

  !!
  !! Test ReLU activation variant
  !!
@Test
  subroutine testReLUActivation()
    type(trainedMLP)            :: mlp_relu
    real(defReal), dimension(3) :: bboxMin, bboxMax, point
    real(defReal)               :: sdf
    real(defReal), parameter    :: TOL = 1.0e-12_defReal

    bboxMin = [-ONE, -ONE, -ONE]
    bboxMax = [ ONE,  ONE,  ONE]
    call mlp_relu % init(3, 2, 2, ACTIVATION_RELU, 0.01_defReal, ONE, bboxMin, bboxMax)
    mlp_relu % weights = mlp % weights
    mlp_relu % biases  = mlp % biases

    ! For input [-0.5, -0.5, -0.5]:
    ! Layer 1: W*[-0.5,-0.5,-0.5] = [-0.5, -0.5]
    ! After ReLU: [0, 0]
    ! Layer 2: W*[0,0] = [0], sdf = tanh(0) = 0
    point = [-0.5_defReal, -0.5_defReal, -0.5_defReal]
    sdf   = mlp_relu % evaluate(point)
    @assertEqual(ZERO, sdf, TOL)

    call mlp_relu % kill()

  end subroutine testReLUActivation

end module trainedMLP_test
