!!
!! Standalone test driver for mlpInference_mod and mlpWeightIO_mod
!!
!! Exercises the trainedMLP forward pass and binary weight file round-trip.
!! Run with: ./Build/mlpInference_sandbox.out
!!
!! This is a temporary test driver. When pFUnit is available, the proper
!! unit tests are in Geometry/Surfaces/Tests/mlpInference_test.f90.
!!
program mlpInference_sandbox

  use numPrecision
  use mlpInference_mod, only : trainedMLP, ACTIVATION_LEAKYRELU, ACTIVATION_RELU
  use mlpWeightIO_mod,  only : writeMLPWeights, readMLPWeights

  implicit none

  type(trainedMLP) :: mlp
  real(defReal), dimension(3) :: bboxMin, bboxMax, point
  real(defReal) :: sdf, expected
  integer :: nPass, nFail
  real(defReal), parameter :: TOL = 1.0e-12_defReal

  nPass = 0
  nFail = 0

  ! ----- Setup: 2-layer MLP (inputDim=3, hiddenDim=2, numLayers=2) -----
  ! Layer 1: W = [[1,0,0],[0,1,0]], b = [0,0]  (projects first 2 components)
  ! Layer 2: W = [[1,1]],           b = [0]    (sums hidden neurons)
  ! Activation: LeakyReLU(0.01), sdfScale = 1.0, bbox = [-1,1]^3
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

  mlp % weights(1, 1, 1) = ONE;  mlp % weights(1, 2, 1) = ZERO; mlp % weights(1, 3, 1) = ZERO
  mlp % weights(2, 1, 1) = ZERO; mlp % weights(2, 2, 1) = ONE;  mlp % weights(2, 3, 1) = ZERO
  mlp % biases(1:2, 1) = ZERO
  mlp % weights(1, 1, 2) = ONE; mlp % weights(1, 2, 2) = ONE
  mlp % biases(1, 2) = ZERO

  print *, '=== mlpInference_mod sandbox tests ==='
  print *

  ! ---- Test 1: isInit flag ----
  call check('isInit set after init()', mlp % isInit .eqv. .true., nPass, nFail)
  call check('inputDim = 3', mlp % inputDim == 3, nPass, nFail)
  call check('weights allocated', allocated(mlp % weights), nPass, nFail)

  ! ---- Test 2: positive forward pass ----
  ! input [0.5, 0.5, 0.5], bbox=[-1,1] so x_norm=x
  ! Layer 1: [0.5, 0.5], LeakyReLU -> [0.5, 0.5]
  ! Layer 2: [1.0], sdf = tanh(1.0) * 1.0
  point    = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  sdf      = mlp % evaluate(point)
  expected = tanh(ONE)
  call checkReal('forward pass (positive inputs)', sdf, expected, TOL, nPass, nFail)

  ! ---- Test 3: normalisation — bbox centre maps to zero output ----
  ! x_norm = 2*(0-(-1))/(1-(-1))-1 = 0
  ! Layer 1: [0, 0] -> LeakyReLU -> [0, 0]
  ! Layer 2: [0], sdf = tanh(0) = 0
  point    = [ZERO, ZERO, ZERO]
  sdf      = mlp % evaluate(point)
  expected = ZERO
  call checkReal('normalisation (centre of bbox -> sdf=0)', sdf, expected, TOL, nPass, nFail)

  ! ---- Test 4: LeakyReLU negative slope ----
  ! input [-0.5, -0.5, -0.5] -> x_norm = [-0.5, -0.5, -0.5]
  ! Layer 1: [-0.5, -0.5], LeakyReLU -> [-0.005, -0.005]
  ! Layer 2: [-0.01], sdf = tanh(-0.01)
  point    = [-0.5_defReal, -0.5_defReal, -0.5_defReal]
  sdf      = mlp % evaluate(point)
  expected = tanh(-0.01_defReal)
  call checkReal('LeakyReLU negative slope', sdf, expected, TOL, nPass, nFail)

  ! ---- Test 5: tanh output bounded in [-1, 1] for extreme inputs ----
  ! tanh saturates to exactly 1.0 or -1.0 in floating point, so use <= not <
  point = [1000.0_defReal, 1000.0_defReal, 1000.0_defReal]
  sdf   = mlp % evaluate(point)
  call check('tanh output <= 1.0 (saturated positive)', sdf <= ONE, nPass, nFail)
  point = [-1000.0_defReal, -1000.0_defReal, -1000.0_defReal]
  sdf   = mlp % evaluate(point)
  call check('tanh output >= -1.0 (saturated negative)', sdf >= -ONE, nPass, nFail)

  ! ---- Test 6: sdfScale multiplier ----
  ! Reinit with sdfScale=5.0, same weights
  ! sdf_scaled = tanh(1.0) * 5.0 for point [0.5, 0.5, 0.5]
  call mlp % kill()
  call mlp % init(3, 2, 2, ACTIVATION_LEAKYRELU, 0.01_defReal, 5.0_defReal, bboxMin, bboxMax)
  mlp % weights(1, 1, 1) = ONE;  mlp % weights(1, 2, 1) = ZERO; mlp % weights(1, 3, 1) = ZERO
  mlp % weights(2, 1, 1) = ZERO; mlp % weights(2, 2, 1) = ONE;  mlp % weights(2, 3, 1) = ZERO
  mlp % biases(1:2, 1) = ZERO
  mlp % weights(1, 1, 2) = ONE; mlp % weights(1, 2, 2) = ONE
  mlp % biases(1, 2) = ZERO

  point    = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  sdf      = mlp % evaluate(point)
  expected = tanh(ONE) * 5.0_defReal
  call checkReal('sdfScale=5.0 scales output', sdf, expected, TOL, nPass, nFail)

  ! ---- Test 7: ReLU variant (negative input -> 0 output) ----
  call mlp % kill()
  call mlp % init(3, 2, 2, ACTIVATION_RELU, 0.01_defReal, ONE, bboxMin, bboxMax)
  mlp % weights(1, 1, 1) = ONE;  mlp % weights(1, 2, 1) = ZERO; mlp % weights(1, 3, 1) = ZERO
  mlp % weights(2, 1, 1) = ZERO; mlp % weights(2, 2, 1) = ONE;  mlp % weights(2, 3, 1) = ZERO
  mlp % biases(1:2, 1) = ZERO
  mlp % weights(1, 1, 2) = ONE; mlp % weights(1, 2, 2) = ONE
  mlp % biases(1, 2) = ZERO

  ! input [-0.5,-0.5,-0.5]: Layer 1 -> [-0.5,-0.5], ReLU -> [0,0], Layer 2 -> [0]
  point    = [-0.5_defReal, -0.5_defReal, -0.5_defReal]
  sdf      = mlp % evaluate(point)
  expected = ZERO
  call checkReal('ReLU: negative input -> zero output', sdf, expected, TOL, nPass, nFail)

  ! ---- Test 8: kill() resets state ----
  call mlp % kill()
  call check('isInit false after kill()', mlp % isInit .eqv. .false., nPass, nFail)
  call check('weights deallocated after kill()', .not. allocated(mlp % weights), nPass, nFail)
  call check('inputDim=0 after kill()', mlp % inputDim == 0, nPass, nFail)

  ! ======================================================
  ! Round-trip tests: write binary file, read it back
  ! ======================================================
  print *
  print *, '=== mlpWeightIO_mod round-trip tests ==='
  print *

  ! Setup a fresh LeakyReLU MLP for round-trip
  call mlp % kill()
  call mlp % init(3, 2, 2, ACTIVATION_LEAKYRELU, 0.01_defReal, ONE, bboxMin, bboxMax)
  mlp % weights(1, 1, 1) = ONE;  mlp % weights(1, 2, 1) = ZERO; mlp % weights(1, 3, 1) = ZERO
  mlp % weights(2, 1, 1) = ZERO; mlp % weights(2, 2, 1) = ONE;  mlp % weights(2, 3, 1) = ZERO
  mlp % biases(1:2, 1) = ZERO
  mlp % weights(1, 1, 2) = ONE; mlp % weights(1, 2, 2) = ONE
  mlp % biases(1, 2) = ZERO

  ! Choose a test point and compute expected output before writing
  point    = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  expected = mlp % evaluate(point)   ! = tanh(1.0)

  ! ---- Test: write binary file ----
  call writeMLPWeights(mlp, '/tmp/test_mlp.bin', point, expected)
  call check('writeMLPWeights does not crash', .true., nPass, nFail)

  ! ---- Test: read binary file back ----
  block
    type(trainedMLP) :: mlp2
    real(defReal)    :: sdf2
    call readMLPWeights(mlp2, '/tmp/test_mlp.bin')
    call check('readMLPWeights: isInit after read', mlp2 % isInit, nPass, nFail)
    call check('readMLPWeights: inputDim', mlp2 % inputDim == 3, nPass, nFail)
    call check('readMLPWeights: hiddenDim', mlp2 % hiddenDim == 2, nPass, nFail)
    call check('readMLPWeights: numLayers', mlp2 % numLayers == 2, nPass, nFail)
    ! Evaluate at same point — should match exactly since same weights
    sdf2 = mlp2 % evaluate(point)
    call checkReal('round-trip: evaluate matches original', sdf2, expected, TOL, nPass, nFail)
    ! Weights preserved
    call checkReal('round-trip: W(1,1,1)', mlp2 % weights(1,1,1), ONE,  TOL, nPass, nFail)
    call checkReal('round-trip: W(1,2,1)', mlp2 % weights(1,2,1), ZERO, TOL, nPass, nFail)
    call checkReal('round-trip: W(1,1,2)', mlp2 % weights(1,1,2), ONE,  TOL, nPass, nFail)
    ! Biases preserved
    call checkReal('round-trip: bias(1,1)', mlp2 % biases(1,1), ZERO, TOL, nPass, nFail)
    call mlp2 % kill()
  end block

  ! ---- Test: evaluate different point after round-trip ----
  block
    type(trainedMLP) :: mlp3
    real(defReal)    :: sdf_orig, sdf_loaded
    real(defReal), dimension(3) :: pt2
    pt2 = [-0.3_defReal, 0.7_defReal, 0.1_defReal]
    sdf_orig   = mlp % evaluate(pt2)
    call readMLPWeights(mlp3, '/tmp/test_mlp.bin')
    sdf_loaded = mlp3 % evaluate(pt2)
    call checkReal('round-trip: evaluate at second point', sdf_loaded, sdf_orig, TOL, nPass, nFail)
    call mlp3 % kill()
  end block

  call mlp % kill()

  ! ---- Summary ----
  print *
  print '(A,I0,A,I0)', '  Passed: ', nPass, '  Failed: ', nFail
  if (nFail == 0) then
    print *, '  ALL TESTS PASSED'
    stop 0
  else
    print *, '  SOME TESTS FAILED'
    stop 1
  end if

contains

  subroutine check(name, condition, nPass, nFail)
    character(*), intent(in)    :: name
    logical, intent(in)         :: condition
    integer, intent(inout)      :: nPass, nFail
    if (condition) then
      print '(A,A)', '  [PASS] ', name
      nPass = nPass + 1
    else
      print '(A,A)', '  [FAIL] ', name
      nFail = nFail + 1
    end if
  end subroutine check

  subroutine checkReal(name, got, expected, tol, nPass, nFail)
    character(*), intent(in)    :: name
    real(defReal), intent(in)   :: got, expected, tol
    integer, intent(inout)      :: nPass, nFail
    if (abs(got - expected) <= tol) then
      print '(A,A)', '  [PASS] ', name
      nPass = nPass + 1
    else
      print '(A,A)', '  [FAIL] ', name
      print '(A,ES20.12,A,ES20.12)', '         got=', got, '  expected=', expected
      nFail = nFail + 1
    end if
  end subroutine checkReal

end program mlpInference_sandbox
