module misclassClerk_class

  use numPrecision
  use tallyCodes
  use display_func,          only : statusMsg
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle
  use outputFile_class,      only : outputFile
  use scoreMemory_class,     only : scoreMemory
  use tallyClerk_inter,      only : tallyClerk, kill_super => kill

  ! Nuclear Data interface
  use nuclearDatabase_inter, only : nuclearDatabase

  ! Material menu -- convert a material name to its matIdx
  use materialMenu_mod,      only : mm_matIdx => matIdx

  ! Reference region: an arbitrary SCONE surface, or logical combination of
  ! surfaces, built with the same surfaceShelf/cell machinery used to define
  ! universes.
  use surfaceShelf_class,    only : surfaceShelf
  use cell_inter,            only : cell
  use cellFactory_func,      only : new_cell_ptr

  implicit none
  private

  !!
  !! Halfspace-misclassification diagnostic for a surrogate (trained) surface
  !!
  !! At every real/virtual collision, compares the material the particle was
  !! actually assigned by the transport geometry (p % matIdx()) against an
  !! independently-evaluated reference region (any SCONE surface, or logical
  !! combination of surfaces, via the usual surfaces+cell machinery). Scores
  !! two counts every report -- total collisions seen, and how many of those
  !! were misclassified -- so their ratio (printed by display(), and left for
  !! post-processing from Res[1]/Res[2] in the output file) is a genuine
  !! rate, independent of whatever population normalisation the physics
  !! package applies (it cancels in the ratio).
  !!
  !! Referencing a surrogate surface's own ground-truth geometry back at
  !! itself should always yield 0% misclassification.
  !!
  !! Private Members:
  !!   insideMatIdx  -> matIdx that represents "inside" the surrogate surface
  !!   flip          -> Invert the reference region's sense (for surrogate
  !!                     models with an inverted halfspace convention)
  !!   handleVirtual -> Score on virtual collisions as well as real ones
  !!   refSurfs      -> Surfaces used to build the reference region
  !!   refCell       -> Cell (any type, any logical combination) defining the
  !!                     reference region from refSurfs
  !!
  !! Interface:
  !!   tallyClerk Interface
  !!
  !! Sample dictionary input:
  !!   misclass { type misclassClerk;
  !!              insideMat fuel;
  !!              # flip 0; handleVirtual 1; #
  !!              surfaces { ref { type sphere; id 1; origin (0 0 0); radius 6.0; } }
  !!              refCell { type simpleCell; id 1; surfaces (-1); }
  !!            }
  !!
  type, public, extends(tallyClerk) :: misclassClerk
    private
    integer(shortInt)    :: insideMatIdx  = 0
    logical(defBool)     :: flip          = .false.
    logical(defBool)     :: handleVirtual = .true.
    type(surfaceShelf)   :: refSurfs
    class(cell), pointer :: refCell => null()
  contains
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportInColl

    ! Output procedures
    procedure :: display
    procedure :: print

  end type misclassClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(misclassClerk), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    character(nameLen), intent(in)      :: name
    character(nameLen)                  :: matName

    ! Assign name
    call self % setName(name)

    ! Material representing "inside" the surrogate surface
    call dict % get(matName, 'insideMat')
    self % insideMatIdx = mm_matIdx(matName)

    ! Optional sign flip and virtual-collision handling
    call dict % getOrDefault(self % flip, 'flip', .false.)
    call dict % getOrDefault(self % handleVirtual, 'handleVirtual', .true.)

    ! Build the reference region from ordinary SCONE surfaces + a cell
    call self % refSurfs % init(dict % getDictPtr('surfaces'))
    self % refCell => new_cell_ptr(dict % getDictPtr('refCell'), self % refSurfs)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(misclassClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    self % insideMatIdx  = 0
    self % flip          = .false.
    self % handleVirtual = .true.

    if (associated(self % refCell)) then
      call self % refCell % kill()
      deallocate(self % refCell)
    end if
    call self % refSurfs % kill()

  end subroutine kill

  !!
  !! Returns array of codes that represent different reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(misclassClerk), intent(in)           :: self
    integer(shortInt), dimension(:), allocatable :: validCodes

    validCodes = [inColl_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(misclassClerk), intent(in) :: self
    integer(shortInt)                :: S

    S = 2

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(misclassClerk), intent(inout)   :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    logical(defBool), intent(in)          :: virtual
    logical(defBool)                      :: hs, ref_hs
    integer(longInt)                      :: addr

    ! Return if collision is virtual but virtual collision handling is off
    if ((.not. self % handleVirtual) .and. virtual) return

    ! Reference region's halfspace at the particle's current position
    ref_hs = .not. self % refCell % inside(p % rGlobal(), p % dirGlobal())
    if (self % flip) ref_hs = .not. ref_hs

    ! Surrogate surface's halfspace, read off the material the transport
    ! geometry already assigned (fresh every collision -- see geom % teleport())
    hs = (p % matIdx() /= self % insideMatIdx)

    ! Bin 1: total collisions seen. Bin 2: how many were misclassified.
    addr = self % getMemAddress()
    call mem % score(ONE, addr)
    if (hs .neqv. ref_hs) call mem % score(ONE, addr + 1_longInt)

  end subroutine reportInColl

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(misclassClerk), intent(in) :: self
    type(scoreMemory), intent(in)    :: mem
    real(defReal) :: totVal, totStd, misVal, misStd, rate
    integer(longInt) :: addr

    addr = self % getMemAddress()
    call mem % getResult(totVal, totStd, addr)
    call mem % getResult(misVal, misStd, addr + 1_longInt)

    rate = ZERO
    if (totVal > ZERO) rate = 100.0_defReal * misVal / totVal

    print '(A,A,A,F8.4,A)', '  ', trim(self % getName()), &
          ' misclassification rate: ', rate, ' %'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(misclassClerk), intent(in) :: self
    class(outputFile), intent(inout) :: outFile
    type(scoreMemory), intent(in)    :: mem
    real(defReal)                    :: val, std
    character(nameLen)               :: name
    integer(longInt)                 :: addr

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Write results: Res(1) = total collisions seen, Res(2) = misclassified.
    ! Their ratio is the misclassification rate.
    addr = self % getMemAddress()
    name = 'Res'
    call outFile % startArray(name, [2])
    call mem % getResult(val, std, addr)
    call outFile % addResult(val, std)
    call mem % getResult(val, std, addr + 1_longInt)
    call outFile % addResult(val, std)
    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

end module misclassClerk_class
