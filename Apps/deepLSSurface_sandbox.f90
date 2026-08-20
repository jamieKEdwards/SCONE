!!
!! Standalone test driver for deepLSSurface_class (paper-faithful DeepLS:
!! shared decoder + per-voxel latent codes)
!!
!! Prerequisites:
!!   python3 scripts/neuralSurface/train_deepls_shared.py \
!!       --input data/rbsphere_unit_train.bin --output /tmp/deepls_shared_smoke.bin \
!!       --nvox 4 4 4 --latent-dim 8 --hidden-dim 16 --num-layers 2 --epochs 20 --device cpu
!!
!! Run:
!!   cmake --build Build --target deepLSSurface_sandbox.out
!!   ./Build/deepLSSurface_sandbox.out
!!
program deepLSSurface_sandbox

  use numPrecision
  use universalVariables, only : INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use deepLSSurface_class, only : deepLSSurface

  implicit none

  type(deepLSSurface)   :: ds
  type(dictionary)            :: dict
  real(defReal)                :: c, d
  real(defReal), dimension(3)  :: r, u
  integer :: nPass, nFail

  nPass = 0
  nFail = 0

  print *, '=== deepLSSurface_class sandbox tests ==='
  print *

  call charToDict(dict, "id 1; weightFile /tmp/deepls_shared_smoke.bin;")
  call ds % init(dict)

  call check('myType = deepLSSurface', ds % myType() == 'deepLSSurface', nPass, nFail)

  block
    real(defReal), dimension(6) :: bb
    bb = ds % boundingBox()
    call checkReal('bbox xmin = -1.5', bb(1), -1.5_defReal, nPass, nFail)
    call checkReal('bbox xmax =  1.5', bb(4),  1.5_defReal, nPass, nFail)
  end block

  r = [0.5_defReal, 0.5_defReal, 0.5_defReal]
  u = [ONE, ZERO, ZERO]
  d = ds % distance(r, u)
  call check('distance returns INF', d >= INF, nPass, nFail)

  print *
  print *, 'evaluate() spot checks (finite-value / no-crash):'
  r = [1.4_defReal, 1.4_defReal, 1.4_defReal]
  c = ds % evaluate(r)
  print '(A,3F6.2,A,ES14.6)', '  r=', r, '  c=', c
  call check('far-corner evaluate is finite', abs(c) < 1.0e10_defReal, nPass, nFail)

  r = [ZERO, ZERO, ZERO]
  c = ds % evaluate(r)
  print '(A,3F6.2,A,ES14.6)', '  r=', r, '  c=', c
  call check('centre evaluate is finite', abs(c) < 1.0e10_defReal, nPass, nFail)

  r = [10.0_defReal, 10.0_defReal, 10.0_defReal]
  c = ds % evaluate(r)
  call check('far outside grid -> positive (outside)', c > ZERO, nPass, nFail)

  call dict % kill()

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

  subroutine checkReal(name, got, expected, nPass, nFail)
    character(*), intent(in)  :: name
    real(defReal), intent(in) :: got, expected
    integer, intent(inout)    :: nPass, nFail
    real(defReal), parameter :: TOL = 1.0e-10_defReal
    if (abs(got - expected) <= TOL) then
      print '(A,A)', '  [PASS] ', name
      nPass = nPass + 1
    else
      print '(A,A)', '  [FAIL] ', name
      print '(A,ES20.12,A,ES20.12)', '         got=', got, '  expected=', expected
      nFail = nFail + 1
    end if
  end subroutine checkReal

end program deepLSSurface_sandbox
