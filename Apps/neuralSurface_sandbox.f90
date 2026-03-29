!!
!! Standalone test driver for neuralSurface_class
!!
!! Tests that neuralSurface correctly loads an MLP from the binary file
!! written by scripts/neuralSurface/roundtrip_test.py and evaluates
!! halfspace membership correctly for a sphere-like MLP.
!!
!! Prerequisites:
!!   python3 scripts/neuralSurface/roundtrip_test.py   (writes /tmp/roundtrip_test.bin)
!!
!! Run:
!!   cmake --build Build --target neuralSurface_sandbox.out
!!   ./Build/neuralSurface_sandbox.out
!!
program neuralSurface_sandbox

  use numPrecision
  use universalVariables, only : INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use neuralSurface_class, only : neuralSurface

  implicit none

  type(neuralSurface)       :: ns
  type(dictionary)          :: dict
  real(defReal)             :: c, d
  real(defReal), dimension(3) :: r, u
  logical(defBool)          :: hs
  integer :: nPass, nFail
  real(defReal), parameter :: TOL = 1.0e-12_defReal

  nPass = 0
  nFail = 0

  print *, '=== neuralSurface_class sandbox tests ==='
  print *

  ! ---- Setup: load MLP from round-trip test binary ----
  ! The binary encodes: W1=[[1,0,0],[0,1,0]], b1=[0,0]
  !                     W2=[[1,1]],           b2=[0]
  ! Activation: LeakyReLU(0.01), sdfScale=1.0, bbox=[-1,1]^3
  !
  ! F(r) = tanh(LeakyReLU(r1) + LeakyReLU(r2)) * 1.0
  !  where r1, r2 are the first two bbox-normalised coordinates.

  call charToDict(dict, "id 1; weightFile /tmp/roundtrip_test.bin;")
  call ns % init(dict)

  ! ---- Test 1: myType ----
  call check('myType = neuralSurface', ns % myType() == 'neuralSurface', nPass, nFail)

  ! ---- Test 2: evaluate at [0.5, 0.5, 0.5] ----
  ! x_norm = [0.5, 0.5, 0.5], LReLU([0.5,0.5]) = [0.5,0.5]
  ! Layer2: 1.0, sdf = tanh(1.0)
  r = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  c = ns % evaluate(r)
  call checkReal('evaluate [0.5,0.5,0.5] = tanh(1.0)', c, tanh(ONE), TOL, nPass, nFail)

  ! ---- Test 3: evaluate at bbox centre [0,0,0] ----
  ! x_norm = [0,0,0] (bbox=[-1,1]), layer1 -> [0,0], LReLU -> [0,0]
  ! Layer2: 0.0, sdf = tanh(0.0) = 0.0
  r = [ZERO, ZERO, ZERO]
  c = ns % evaluate(r)
  call checkReal('evaluate at bbox centre = 0', c, ZERO, TOL, nPass, nFail)

  ! ---- Test 4: halfspace positive (outside) ----
  ! c = tanh(1.0) > 0 => +ve halfspace
  r = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  u = [ONE,         ZERO,        ZERO        ]
  hs = ns % halfspace(r, u)
  call check('halfspace +ve for positive c', hs .eqv. .true., nPass, nFail)

  ! ---- Test 5: halfspace negative (inside) ----
  ! r = [-0.5,-0.5,-0.5]: c = tanh(-0.01) < 0 => -ve halfspace
  r = [-0.5_defReal, -0.5_defReal, -0.5_defReal]
  u = [ ONE,          ZERO,         ZERO        ]
  hs = ns % halfspace(r, u)
  call check('halfspace -ve for negative c', hs .eqv. .false., nPass, nFail)

  ! ---- Test 6: distance stub returns INF ----
  r = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  u = [ONE,         ZERO,        ZERO        ]
  d = ns % distance(r, u)
  call check('distance returns INF', d >= INF, nPass, nFail)

  ! ---- Test 7: boundingBox matches weight file bbox ----
  block
    real(defReal), dimension(6) :: bb
    bb = ns % boundingBox()
    call checkReal('bbox xmin = -1', bb(1), -ONE, TOL, nPass, nFail)
    call checkReal('bbox xmax =  1', bb(4),  ONE, TOL, nPass, nFail)
  end block

  ! ---- Test 8: going() FD direction check ----
  ! At centre (c=0), going in direction u=[1,0,0]:
  ! r + eps * u = [eps, 0, 0], evaluate -> tanh(LeakyReLU(eps)) > 0 => going +ve
  r = [ZERO, ZERO, ZERO]
  u = [ONE,  ZERO, ZERO]
  hs = ns % going(r, u)
  call check('going +ve: direction toward positive halfspace', hs .eqv. .true., nPass, nFail)

  ! ---- Test 9: re-init after kill (init cleans up old MLP) ----
  call dict % kill()
  call charToDict(dict, "id 2; weightFile /tmp/roundtrip_test.bin;")
  call ns % kill()
  call ns % init(dict)
  call check('re-init: myType still neuralSurface', ns % myType() == 'neuralSurface', nPass, nFail)
  r = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  c = ns % evaluate(r)
  call checkReal('re-init: evaluate correct after re-load', c, tanh(ONE), TOL, nPass, nFail)

  call ns % kill()
  call dict % kill()

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

end program neuralSurface_sandbox
