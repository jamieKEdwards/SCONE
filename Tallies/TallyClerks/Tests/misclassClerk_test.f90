module misclassClerk_test

  use numPrecision
  use funit
  use dictionary_class,          only : dictionary
  use dictParser_func,           only : charToDict
  use particle_class,            only : particle
  use scoreMemory_class,         only : scoreMemory
  use outputFile_class,          only : outputFile
  use testNeutronDatabase_class, only : testNeutronDatabase
  use materialMenu_mod,          only : mm_init => init, mm_kill => kill, mm_matIdx => matIdx
  use misclassClerk_class,       only : misclassClerk

  implicit none

  ! Two materials for the material menu: fuel (the surrogate's "inside" material)
  ! and water (everything else)
  character(*), parameter :: MAT_DEF = &
    " fuel  { temp 1; composition {} } &
    & water { temp 1; composition {} } "

  ! Reference region: unit sphere at the origin, no flip
  character(*), parameter :: CLERK_DEF = &
    " type misclassClerk; insideMat fuel; &
    & surfaces { ref { type sphere; id 1; origin (0.0 0.0 0.0); radius 1.0; } } &
    & refCell { type simpleCell; id 1; surfaces (-1); } "

  ! Same reference region, with the sense flipped
  character(*), parameter :: CLERK_DEF_FLIP = &
    " type misclassClerk; insideMat fuel; flip 1; &
    & surfaces { ref { type sphere; id 1; origin (0.0 0.0 0.0); radius 1.0; } } &
    & refCell { type simpleCell; id 1; surfaces (-1); } "

  ! Reference region with virtual collisions ignored
  character(*), parameter :: CLERK_DEF_NOVIRT = &
    " type misclassClerk; insideMat fuel; handleVirtual 0; &
    & surfaces { ref { type sphere; id 1; origin (0.0 0.0 0.0); radius 1.0; } } &
    & refCell { type simpleCell; id 1; surfaces (-1); } "

contains

  !!
  !! Build a particle at position r with a given matIdx
  !!
  function buildParticle(r, matIdx) result(p)
    real(defReal), dimension(3), intent(in) :: r
    integer(shortInt), intent(in)           :: matIdx
    type(particle)                          :: p

    call p % build(r, [ONE, ZERO, ZERO], ONE, ONE)
    call p % setMatIdx(matIdx)

  end function buildParticle

  !!
  !! Test scoring of matched and mismatched cases, inside and outside the
  !! reference region
  !!
@Test
  subroutine testScoring()
    type(dictionary)           :: matDict, clerkDict
    type(misclassClerk)        :: clerk
    type(scoreMemory)          :: mem
    type(particle)              :: p
    type(testNeutronDatabase)  :: nucData
    type(outputFile)            :: outF
    integer(shortInt)          :: fuelIdx, waterIdx
    real(defReal)               :: total, totalSTD, mis, misSTD
    real(defReal), parameter    :: TOL = 1.0E-9

    call charToDict(matDict, MAT_DEF)
    call mm_init(matDict)
    fuelIdx  = mm_matIdx('fuel')
    waterIdx = mm_matIdx('water')

    call charToDict(clerkDict, CLERK_DEF)
    call clerk % init(clerkDict, 'misclass')

    call mem % init(int(clerk % getSize(), longInt), 1)
    call clerk % setMemAddress(1_longInt)

    call nucData % build(0.3_defReal)

    ! Matched: inside the sphere, assigned the "inside" material
    p = buildParticle([ZERO, ZERO, ZERO], fuelIdx)
    call clerk % reportInColl(p, nucData, mem, .false.)

    ! Matched: outside the sphere, assigned the "outside" material
    p = buildParticle([2.0_defReal, ZERO, ZERO], waterIdx)
    call clerk % reportInColl(p, nucData, mem, .false.)

    ! Mismatched: inside the sphere, but assigned the "outside" material
    p = buildParticle([ZERO, ZERO, ZERO], waterIdx)
    call clerk % reportInColl(p, nucData, mem, .false.)

    ! Mismatched: outside the sphere, but assigned the "inside" material
    p = buildParticle([2.0_defReal, ZERO, ZERO], fuelIdx)
    call clerk % reportInColl(p, nucData, mem, .false.)

    call mem % reduceBins()
    call mem % closeCycle(ONE)

    call mem % getResult(total, totalSTD, clerk % getMemAddress())
    call mem % getResult(mis, misSTD, clerk % getMemAddress() + 1_longInt)

    @assertEqual(4.0_defReal, total, TOL, 'Total collisions scored')
    @assertEqual(2.0_defReal, mis,   TOL, 'Misclassified count')
    @assertEqual(2, clerk % getSize(), 'Memory size')

    ! Verify output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call clerk % print(outF, mem)
    @assertTrue(outF % isValid())

    call clerk % kill()
    call nucData % kill()
    call clerkDict % kill()
    call matDict % kill()
    call mm_kill()

  end subroutine testScoring

  !!
  !! Test that "flip" inverts the sense of the reference region
  !!
@Test
  subroutine testFlip()
    type(dictionary)           :: matDict, clerkDict
    type(misclassClerk)        :: clerk
    type(scoreMemory)          :: mem
    type(particle)              :: p
    type(testNeutronDatabase)  :: nucData
    integer(shortInt)          :: fuelIdx, waterIdx
    real(defReal)               :: total, totalSTD, mis, misSTD
    real(defReal), parameter    :: TOL = 1.0E-9

    call charToDict(matDict, MAT_DEF)
    call mm_init(matDict)
    fuelIdx  = mm_matIdx('fuel')
    waterIdx = mm_matIdx('water')

    call charToDict(clerkDict, CLERK_DEF_FLIP)
    call clerk % init(clerkDict, 'misclass')

    call mem % init(int(clerk % getSize(), longInt), 1)
    call clerk % setMemAddress(1_longInt)

    call nucData % build(0.3_defReal)

    ! Both of these are correctly matched WITHOUT a flip -- with flip they
    ! must both come out as misclassified
    p = buildParticle([ZERO, ZERO, ZERO], fuelIdx)
    call clerk % reportInColl(p, nucData, mem, .false.)

    p = buildParticle([2.0_defReal, ZERO, ZERO], waterIdx)
    call clerk % reportInColl(p, nucData, mem, .false.)

    call mem % reduceBins()
    call mem % closeCycle(ONE)

    call mem % getResult(total, totalSTD, clerk % getMemAddress())
    call mem % getResult(mis, misSTD, clerk % getMemAddress() + 1_longInt)

    @assertEqual(2.0_defReal, total, TOL, 'Total collisions scored')
    @assertEqual(2.0_defReal, mis,   TOL, 'Misclassified count with flip')

    call clerk % kill()
    call nucData % kill()
    call clerkDict % kill()
    call matDict % kill()
    call mm_kill()

  end subroutine testFlip

  !!
  !! Test that handleVirtual = 0 causes virtual collisions to be ignored
  !!
@Test
  subroutine testHandleVirtual()
    type(dictionary)           :: matDict, clerkDict
    type(misclassClerk)        :: clerk
    type(scoreMemory)          :: mem
    type(particle)              :: p
    type(testNeutronDatabase)  :: nucData
    integer(shortInt)          :: fuelIdx
    real(defReal)               :: total, totalSTD
    real(defReal), parameter    :: TOL = 1.0E-9

    call charToDict(matDict, MAT_DEF)
    call mm_init(matDict)
    fuelIdx = mm_matIdx('fuel')

    call charToDict(clerkDict, CLERK_DEF_NOVIRT)
    call clerk % init(clerkDict, 'misclass')

    call mem % init(int(clerk % getSize(), longInt), 1)
    call clerk % setMemAddress(1_longInt)

    call nucData % build(0.3_defReal)

    ! Virtual collision -- must be ignored
    p = buildParticle([ZERO, ZERO, ZERO], fuelIdx)
    call clerk % reportInColl(p, nucData, mem, .true.)

    ! Real collision -- must be scored
    p = buildParticle([ZERO, ZERO, ZERO], fuelIdx)
    call clerk % reportInColl(p, nucData, mem, .false.)

    call mem % reduceBins()
    call mem % closeCycle(ONE)

    call mem % getResult(total, totalSTD, clerk % getMemAddress())

    @assertEqual(1.0_defReal, total, TOL, 'Virtual collision should not be scored')

    call clerk % kill()
    call nucData % kill()
    call clerkDict % kill()
    call matDict % kill()
    call mm_kill()

  end subroutine testHandleVirtual

end module misclassClerk_test
