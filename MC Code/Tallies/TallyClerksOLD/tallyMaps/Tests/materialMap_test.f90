module materialMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particle
  use dictionary_class,        only : dictionary
  use nuclearDataRegistry_mod, only : build_NuclearData, kill_nuclearData

  use materialMap_class,       only : materialMap

  implicit none


@testCase
  type, extends(TestCase) :: test_materialMap
    private
    type(materialMap),allocatable :: map_noUndef
    type(materialMap),allocatable :: map_Undef
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_materialMap


  !!
  !! Test parameters
  !!
  character(*),dimension(*),parameter :: MAT_NAMES=['mat1','mat2','mat3','mat4','mat5']
  character(*),dimension(*),parameter :: MAT_IN_MAP =['mat2','mat3','mat5']




contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_materialMap), intent(inout) :: this
    integer(shortInt)                 :: temp
    integer(shortInt)                 :: i
    type(dictionary)                  :: dict, tempDict1, tempDict2, tempDict3
    type(dictionary)                  :: mapDict1

    ! Initialise dictionaries
    call dict % init(6)
    call tempDict1 % init(6)
    call tempDict2 % init(6)
    call tempDict3 % init(6)
    call mapDict1 % init(2)

    !*** Provisional -> allow registry without handles
    call tempDict2 % store('myMat','datalessMaterials')

    ! Store empty- handles dictionary
    call dict % store('handles', tempDict2)

    ! Create materials dictionary of empty dictionaries
    do i=1,size(MAT_NAMES)
      call tempDict1 % store(MAT_NAMES(i), tempDict3)

    end do

    ! Store materials dictionary in nuclearData dict
    call dict % store('materials',tempDict1)

    ! Build nuclear data
    call build_NuclearData(dict)

    ! Build material map definition
    call mapDict1 % store('materials', MAT_IN_MAP)
    allocate(this % map_noUndef, source = materialMap(mapDict1))

    call mapDict1 % store('undefBin','true')
    allocate(this % map_undef, source = materialMap(mapDict1))

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_materialMap), intent(inout) :: this

    call kill_nuclearData()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Mapping test without undefined bin
  !!
@Test
  subroutine testMappingNoUndefined(this)
    class(test_materialMap), intent(inout)   :: this
    type(particle)                           :: p
    integer(shortInt)                        :: i
    integer(shortInt),dimension(5)           :: bins
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [0, 1, 2, 0, 3]

    do i=1,5
      call p % setMatIdx(i)
      bins(i) = this % map_noUndef % map(p)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMappingNoUndefined


  !!
  !! Mapping test with undefined bin
  !!
@Test
  subroutine testMappingUndefined(this)
    class(test_materialMap), intent(inout)   :: this
    type(particle)                           :: p
    integer(shortInt)                        :: i
    integer(shortInt),dimension(5)           :: bins
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [4, 1, 2, 4, 3]

    do i=1,5
      call p % setMatIdx(i)
      bins(i) = this % map_undef % map(p)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMappingUndefined


  !!
  !! Test number of bins inquiry
  !!
@Test
  subroutine testNumberOfBinsInquiry(this)
    class(test_materialMap), intent(inout) :: this

    @assertEqual(3, this % map_noUndef % bins(),'materialMap without undefined bin')
    @assertEqual(4, this % map_undef % bins(),'materialMap with undefined bin')

  end subroutine testNumberOfBinsInquiry


end module materialMap_test