!!
!! Python → Fortran round-trip test
!!
!! Reads the binary weight file written by scripts/neuralSurface/roundtrip_test.py
!! and verifies the loaded MLP gives the expected output.
!!
!! Run:
!!   python3 scripts/neuralSurface/roundtrip_test.py   (writes /tmp/roundtrip_test.bin)
!!   cmake --build Build --target mlpRoundtrip_sandbox.out
!!   ./Build/mlpRoundtrip_sandbox.out
!!
program mlpRoundtrip_sandbox

  use numPrecision
  use mlpInference_mod, only : trainedMLP, ACTIVATION_LEAKYRELU
  use mlpWeightIO_mod,  only : readMLPWeights

  implicit none

  type(trainedMLP) :: mlp
  real(defReal)    :: sdf
  real(defReal), dimension(3) :: point
  integer :: nPass, nFail
  real(defReal), parameter :: TOL = 1.0e-12_defReal

  nPass = 0
  nFail = 0

  print *, '=== Python -> Fortran round-trip test ==='
  print *

  ! Read the file written by roundtrip_test.py
  call readMLPWeights(mlp, '/tmp/roundtrip_test.bin')

  ! ---- Architecture checks ----
  call check('isInit after readMLPWeights', mlp % isInit,          nPass, nFail)
  call check('inputDim = 3',               mlp % inputDim  == 3,  nPass, nFail)
  call check('hiddenDim = 2',              mlp % hiddenDim == 2,  nPass, nFail)
  call check('numLayers = 2',              mlp % numLayers == 2,  nPass, nFail)
  call check('activation = LeakyReLU',                            &
             mlp % activationType == ACTIVATION_LEAKYRELU, nPass, nFail)

  ! ---- Weight checks ----
  ! Layer 1 W = [[1,0,0],[0,1,0]]: W(out,in,layer)
  call checkReal('W(1,1,1) = 1.0', mlp % weights(1,1,1), ONE,  TOL, nPass, nFail)
  call checkReal('W(1,2,1) = 0.0', mlp % weights(1,2,1), ZERO, TOL, nPass, nFail)
  call checkReal('W(1,3,1) = 0.0', mlp % weights(1,3,1), ZERO, TOL, nPass, nFail)
  call checkReal('W(2,1,1) = 0.0', mlp % weights(2,1,1), ZERO, TOL, nPass, nFail)
  call checkReal('W(2,2,1) = 1.0', mlp % weights(2,2,1), ONE,  TOL, nPass, nFail)
  call checkReal('W(2,3,1) = 0.0', mlp % weights(2,3,1), ZERO, TOL, nPass, nFail)
  ! Layer 2 W = [[1,1]]
  call checkReal('W(1,1,2) = 1.0', mlp % weights(1,1,2), ONE,  TOL, nPass, nFail)
  call checkReal('W(1,2,2) = 1.0', mlp % weights(1,2,2), ONE,  TOL, nPass, nFail)

  ! ---- Forward pass: [0.5, 0.5, 0.5] should give tanh(1.0) ----
  ! Layer 1: x_norm = [0.5,0.5,0.5], h = LeakyReLU([0.5,0.5]) = [0.5,0.5]
  ! Layer 2: z = 1.0, sdf = tanh(1.0) * 1.0
  point = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  sdf   = mlp % evaluate(point)
  call checkReal('evaluate([0.5,0.5,0.5]) = tanh(1.0)', sdf, tanh(ONE), TOL, nPass, nFail)

  ! ---- Second point: [-0.5, -0.5, -0.5] should give tanh(-0.01) ----
  ! LeakyReLU([-0.5,-0.5]) = [-0.005,-0.005], z = -0.01, sdf = tanh(-0.01)
  point = [-0.5_defReal, -0.5_defReal, -0.5_defReal]
  sdf   = mlp % evaluate(point)
  call checkReal('evaluate([-0.5,-0.5,-0.5]) = tanh(-0.01)', sdf, tanh(-0.01_defReal), TOL, nPass, nFail)

  call mlp % kill()
  call check('isInit false after kill()', .not. mlp % isInit, nPass, nFail)

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
    character(*), intent(in) :: name
    logical,      intent(in) :: condition
    integer, intent(inout)   :: nPass, nFail
    if (condition) then
      print '(A,A)', '  [PASS] ', name
      nPass = nPass + 1
    else
      print '(A,A)', '  [FAIL] ', name
      nFail = nFail + 1
    end if
  end subroutine check

  subroutine checkReal(name, got, expected, tol, nPass, nFail)
    character(*), intent(in)  :: name
    real(defReal), intent(in) :: got, expected, tol
    integer, intent(inout)    :: nPass, nFail
    if (abs(got - expected) <= tol) then
      print '(A,A)', '  [PASS] ', name
      nPass = nPass + 1
    else
      print '(A,A)', '  [FAIL] ', name
      print '(A,ES20.12,A,ES20.12)', '         got=', got, '  expected=', expected
      nFail = nFail + 1
    end if
  end subroutine checkReal

end program mlpRoundtrip_sandbox
