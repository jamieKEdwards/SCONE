module neutronCEstd_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG
  use intMap_class,                  only : intMap

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON
  use particleDungeon_class,         only : particleDungeon

  ! Abstarct interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data Interfaces
  use nuclearDataReg_mod,            only : ndReg_getNeutronCE => getNeutronCE
  use nuclearDatabase_inter,         only : nuclearDatabase
  use ceNeutronDatabase_inter,       only : ceNeutronDatabase
  use aceNeutronDatabase_class,      only : aceNeutronDatabase
  use aceNeutronDatabase_class,      only : aceNeutronDatabase_CptrCast
  use ceNeutronMaterial_class,       only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use ceNeutronNuclide_inter,        only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,       only : aceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,       only : aceNeutronNuclide
  use materialMenu_mod,              only : materialItem, mm_getMatPtr => getMatPtr

  ! Nuclear reactions
  use reactionHandle_inter,          only : reactionHandle
  use uncorrelatedReactionCE_inter,  only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use neutronScatter_class,          only : neutronScatter, neutronScatter_TptrCast
  use fissionCE_class,               only : fissionCE, fissionCE_TptrCast

  ! Cross-Section Packages
  use neutronXsPackages_class,       only : neutronMicroXSs

  ! Scattering procedures
  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS, targetVelocity_DBRCXS,  &
                                     asymptoticInelasticScatter, targetVelocity_relE

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for CE neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE    -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE    -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 20.0)
  !!  thresh_E -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!              Target movment is sampled if neutron energy E < kT * thresh_E where
  !!              kT is target material temperature in [MeV]. (default = 400.0)
  !!  thresh_A -> Mass threshold for explicit tratment of target nuclide movement [Mn].
  !!              Target movment is sampled if target mass A < thresh_A. (default = 1.0)
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type             neutronCEstd;
  !!   #minEnergy       <real>;#
  !!   #maxEnergy       <real>;#
  !!   #energyThreshold <real>;#
  !!   #massThreshold   <real>;#
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronCEstd
    private
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(ceNeutronDatabase), pointer, public :: xsData => null()
    class(ceNeutronMaterial), pointer, public :: mat    => null()
    class(ceNeutronNuclide),  pointer, public :: nuc    => null()
    class(aceNeutronNuclide), pointer, public :: aceNuc => null()
    class(aceNeutronDatabase), pointer, public:: aceData=> null()
    integer(longInt), public :: count = 0
    integer(longInt), public :: rejCount = 0
    real(defReal), public :: frac = 0




    !! Settings - private
    real(defReal) :: minE
    real(defReal) :: maxE
    real(defReal) :: tresh_E
    real(defReal) :: tresh_A
    real(defReal) :: DBRCeMin
    real(defReal) :: DBRCeMax

  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: elastic
    procedure :: inelastic
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs

    ! Local procedures
    procedure,private :: scatterFromFixed
    procedure,private :: scatterFromMoving
    procedure,private :: scatterInLAB
  end type neutronCEstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronCEstd), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter :: Here = 'init (neutronCEstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Read settings for neutronCEstd
    ! Maximum and minimum energy
    call dict % getOrDefault(self % minE,'minEnergy',1.0E-11_defReal)
    call dict % getOrDefault(self % maxE,'maxEnergy',20.0_defReal)

    ! Thermal scattering kernel thresholds
    call dict % getOrDefault(self % tresh_E, 'energyThreshold', 400.0_defReal)
    call dict % getOrDefault(self % tresh_A, 'massThreshold', 1.0_defReal)

    ! Verify settings
    if( self % minE < ZERO ) call fatalError(Here,'-ve minEnergy')
    if( self % maxE < ZERO ) call fatalError(Here,'-ve maxEnergy')
    if( self % minE >= self % maxE) call fatalError(Here,'minEnergy >= maxEnergy')
    if( self % tresh_E < 0) call fatalError(Here,' -ve energyThreshold')
    if( self % tresh_A < 0) call fatalError(Here,' -ve massThreshold')

    !Read in DBRC max and min energy cut off
    !if dict % isPresent(self % DBRC) then
    ! maybe can go without the if function
    call dict % getOrDefault(self % DBRCeMin,'DBRCeMin', (1.0E-8_defReal))
    call dict % getOrDefault(self % DBRCeMax,'DBRCeMax', (200E-6_defReal))

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMicroXSs)                :: microXSs
    real(defReal)                        :: r, rel_E
    character(100),parameter :: Here = 'sampleCollision (neutronCEstd_class.f90)'

    ! Verify that particle is CE neutron
    if(p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only CE Neutron. Was given MG '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronCE()
    if(.not.associated(self % xsData)) call fatalError(Here, 'There is no active Neutron CE data!')

    ! Verify and load material pointer
    self % mat => ceNeutronMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, 'Material is not ceNeutronMaterial')

    ! Select collision nuclide
    call self % mat % sampleNuclide(p % E, p % pRNG, collDat % nucIdx, rel_E)

    ! If nuclide of index -2 is returned then TMS collision site and nuclide is rejected
    if ((self % mat % matHasTMS) .and. (collDat % nucIdx == -2)) then
      collDat % MT = noInteraction
      return
    end if

    self % nuc => ceNeutronNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
    if(.not.associated(self % mat)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

    ! Find relative energy micro cross sections for sampling reaction
    if (self % mat % matHasTMS) then
      call self % nuc % getMicroXSs(microXss, rel_E, p % pRNG)
    else
      call self % nuc % getMicroXSs(microXss, p % E, p % pRNG)
    end if

    ! Sample reaction
    r = p % pRNG % get()
    collDat % MT = microXss % invert(r)

  end subroutine sampleCollision

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(fissionCE), pointer             :: fission
    type(materialItem), pointer          :: matItem
    type(neutronMicroXSs)                :: microXSs
    type(particleState)                  :: pTemp
    real(defReal),dimension(3)           :: r, dir
    integer(shortInt)                    :: n, i
    real(defReal)                        :: wgt, w0, rand1, E_out, mu, phi
    real(defReal)                        :: sig_nufiss, sig_tot, k_eff
    real(defReal)                        :: kT, matkT, deltakT, rel_E, A
    character(100),parameter             :: Here = 'implicit (neutronCEstd_class.f90)'

    ! Generate fission sites if nuclide is fissile
    if ( self % nuc % isFissile()) then
      ! Obtain required data
      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites


      matItem => mm_getMatPtr(p % matIdx())

      !! If material uses TMS then relative energy needs to be sampled for
      !! the cross section look ups which determine fission sites.
      !!
      !! Althought the information below is to find `broadened' cross sections for TMS, the dopplerCorrectionFactor
      !! is not used as it would cancel out in the sig_nufiss/sig_tot
      !!
      if (matItem % hasTMS) then

        ! Bring through material and nuclide information
        matkT = (kBoltzmann * matItem % T) / joulesPerMeV
        kT = self % nuc % getkT()
        A = self % nuc % getMass()
        deltakT = matkT - kT
        rel_E = targetVelocity_relE(p % E, A, deltakT, p % pRNG)

        ! Avoid sampled relative energy from MB dist extending into energies outside system range
        if (rel_E < self % minE) then
          rel_E = self % minE
        end if

        if (self % maxE < rel_E) then
          rel_E = self % maxE
        end if

        ! Look up cross sections with relative E
        call self % nuc % getMicroXSs(microXSs, rel_E, p % pRNG)
      else
        call self % nuc % getMicroXSs(microXSs, p % E, p % pRNG)
      end if

      sig_nufiss = microXSs % nuFission
      sig_tot    = microXSs % total

      ! Sample number of fission sites generated
      ! Support -ve weight particles
      n = int(abs( (wgt * sig_nufiss) / (w0 * sig_tot * k_eff)) + rand1, shortInt)

      ! Shortcut particle generation if no particles were sampled
      if (n < 1) return

      ! Get fission Reaction
      fission => fissionCE_TptrCast(self % xsData % getReaction(N_FISSION, collDat % nucIdx))
      if(.not.associated(fission)) call fatalError(Here, "Failed to get fissionCE")

      ! Store new sites in the next cycle dungeon
      wgt =  sign(w0, wgt)
      r   = p % rGlobal()

      do i=1,n
        call fission % sampleOut(mu, phi, E_out, p % E, p % pRNG)
        dir = rotateVector(p % dirGlobal(), mu, phi)

        if (E_out > self % maxE) E_out = self % maxE

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        pTemp       = p

        ! Overwrite position, direction, energy and weight
        pTemp % r   = r
        pTemp % dir = dir
        pTemp % E   = E_out
        pTemp % wgt = wgt

        call nextCycle % detain(pTemp)
      end do
    end if

  end subroutine implicit

  !!
  !! Process capture reaction
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine capture

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead = .true.

  end subroutine fission

  !!
!! Process elastic scattering
!!
!! All CE elastic scattering happens in the CM frame
!!
subroutine elastic(self, p, collDat, thisCycle, nextCycle)
  class(neutronCEstd), intent(inout)      :: self
  class(particle), intent(inout)          :: p
  type(collisionData), intent(inout)      :: collDat
  class(particleDungeon),intent(inout)    :: thisCycle
  class(particleDungeon),intent(inout)    :: nextCycle
  class(uncorrelatedReactionCE), pointer  :: reac
  type(materialItem), pointer             :: matItem
  logical(defBool)                        :: isFixed, nucDBRC
  real(defReal)                           :: matkT
  character(100),parameter :: Here = 'elastic (neutronCEstd_class.f90)'

  ! Get reaction
  reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
  if(.not.associated(reac)) call fatalError(Here,'Failed to get elastic neutron scatter')

  ! Get collision data
  collDat % A =  self % nuc % getMass()

  ! Bring through material information
  matItem => mm_getMatPtr(p % matIdx())

  ! Check for TMS
  if (matItem % hasTMS) then

    ! Set collision temp as material TMS temp in MeV
    matkT = (kBoltzmann * matItem % T) / joulesPerMeV
    collDat % kT = matkT

  ! Otherwise use nuc temp
  else
    collDat % kT = self % nuc % getkT()
  end if

  !cast pointer to aceNeutronNuclide
  self % aceNuc => aceNeutronNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
  if(.not.associated(self % aceNuc)) call fatalError(Here, 'Failed to retive ACE Neutron Nuclide')

  !cast pointer to aceNeutronDatabase
  self % aceData => aceNeutronDatabase_CptrCast(self % xsData)
  if(.not.associated(self % aceData)) call fatalError(Here, 'Failed to retive ACE Neutron Database')

  isFixed = (p % E > collDat % kT * self % tresh_E) .and. (collDat % A > self % tresh_A)
  nucDBRC = ( self % aceNuc % isNucDBRC .and. self % aceData % hasDBRC)

  ! Apply criterion for Free-Gas vs Fixed Target scattering
  if (.not. reac % inCMFrame()) then
    call self % scatterInLAB(p, collDat, reac)
  elseif (isFixed .and. .not. nucDBRC) then
    call self % scatterFromFixed(p, collDat, reac)
  else
    call self % scatterFromMoving(p, collDat, reac)
  end if

end subroutine elastic

  !!
  !! Process inelastic scattering
  !!
  subroutine inelastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)     :: self
    class(particle), intent(inout)         :: p
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon),intent(inout)   :: thisCycle
    class(particleDungeon),intent(inout)   :: nextCycle
    class(uncorrelatedReactionCE), pointer :: reac
    character(100),parameter  :: Here =' inelastic (neutronCEstd_class.f90)'

    ! Invert inelastic scattering and Get reaction
    collDat % MT = self % nuc % invertInelastic(p % E, p % pRNG)
    reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if(.not.associated(reac)) call fatalError(Here, "Failed to get scattering reaction")

    ! Scatter particle
    if (reac % inCMFrame()) then
      collDat % A =  self % nuc % getMass()
      call self % scatterFromFixed(p, collDat, reac)
    else
      call self % scatterInLAB(p, collDat, reac)
    end if

    ! Apply weight change
    p % w = p % w * reac % release(p % E)

  end subroutine inelastic

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(neutronCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    if (p % E < self % minE ) p % isDead = .true.

  end subroutine cutoffs

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, p, collDat, reac)
    class(neutronCEstd), intent(inout)        :: self
    class(particle), intent(inout)            :: p
    type(collisionData), intent(inout)        :: collDat
    class(uncorrelatedReactionCE), intent(in) :: reac
    real(defReal)                             :: phi    ! Azimuthal scatter angle
    real(defReal)                             :: E_out, mu

    ! Sample scattering angles and post-collision energy
    call reac % sampleOut(mu, phi, E_out, p % E, p % pRNG)

    ! Update neutron state
    p % E = E_out
    call p % rotate(mu, phi)
    collDat % muL = mu

  end subroutine scatterInLAB

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self, p, collDat, reac)
    class(neutronCEstd), intent(inout)         :: self
    class(particle), intent(inout)             :: p
    type(collisionData), intent(inout)         :: collDat
    class(uncorrelatedReactionCE), intent(in)  :: reac
    real(defReal)                              :: phi
    real(defReal)                              :: E_out
    real(defReal)                              :: E_outCM, mu
    integer(shortInt)                          :: MT

    ! Read data
    MT     = collDat % MT

    ! Sample mu , phi and outgoing energy
    call reac % sampleOut(mu, phi, E_outCM, p % E, p % pRNG)

    ! Save incident energy
    E_out = p % E

    if( MT == N_N_elastic) then
      call asymptoticScatter(E_out, mu, collDat % A)
    else
      call asymptoticInelasticScatter(E_out, mu, E_outCM, collDat % A)
    end if

    ! Update particle state
    call p % rotate(mu, phi)
    p % E = E_out
    collDat % muL = mu

  end subroutine scatterFromFixed

  !!
  !! Subroutine to perform scattering from moving target
  !! Supports only elastic collisions
  !!
  subroutine scatterFromMoving(self, p, collDat, reac)
    class(neutronCEstd), intent(inout)         :: self
    class(particle), intent(inout)             :: p
    type(collisionData),intent(inout)          :: collDat
    class(uncorrelatedReactionCE), pointer     :: reac
    integer(shortInt)                          :: MT, nucIdx
    real(defReal)                              :: A, kT, mu
    real(defReal),dimension(3)                 :: V_n           ! Neutron velocity (vector)
    real(defReal)                              :: U_n           ! Neutron speed (scalar)
    real(defReal),dimension(3)                 :: dir_pre       ! Pre-collision direction
    real(defReal),dimension(3)                 :: dir_post      ! Post-collicion direction
    real(defReal),dimension(3)                 :: V_t, V_cm     ! Target and CM velocity
    real(defReal)                              :: phi, dummy    ! Target and CM velocity
    real(defReal)                              :: TmajXS
    logical(defBool)                           :: eRange, nucDBRC
    character(100), parameter :: Here = 'Scatter From Moving (neutronCEstd_class.f90)'

    nucIdx = 0

    ! Read data
    A      = collDat % A
    kT     = collDat % kT
    MT     = collDat % MT
    nucIdx = collDat % nucIdx

    ! Get neutron direction and velocity
    dir_pre = p % dirGlobal()
    V_n     = dir_pre * sqrt(p % E)


    !cast pointer to aceNeutronNuclide
    self % aceNuc => aceNeutronNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
    if(.not.associated(self % aceNuc)) call fatalError(Here, 'Failed to retive ACE Neutron Nuclide')

    !cast pointer to aceNeutronDatabase
    self % aceData => aceNeutronDatabase_CptrCast(self % xsData)
    if(.not.associated(self % aceData)) call fatalError(Here, 'Failed to retive ACE Neutron Database')



    ! Sample velocity of target with constant xs or with dbrc depending on switch and input files
    ! dbrc requires energy less than emax, greater than emin and nuclide to be listed with 0K XS
    eRange = ((p % E <= self % DBRCeMax) .and. (self % DBRCeMin <= p % E))
    nucDBRC = ( self % aceNuc % isNucDBRC .and. self % aceData % hasDBRC)

    if (eRange .and. nucDBRC) then

      !Retrive 0K nuclide index from dbrc nuclide map
      nucIdx = self % aceData % intMapDBRCnucs % get(nucIdx)

      ! Reassign pointer for the 0K nuclide
      self % aceNuc => aceNeutronNuclide_CptrCast(self % xsData % getNuclide(nucIdx))

      ! set temp majorant
      TmajXS = self % aceData % updateTempMicroMajorantXS(p % E, kT, A, nucIdx)

      ! use dbrc (non constant cross section) to sample target velocity
      V_t = targetVelocity_DBRCXS(self % aceNuc, p % E, dir_pre, A, kT, p % pRNG, TmajXS)

    else
      V_t = targetVelocity_constXS(p % E, dir_pre, A, kT, p % pRNG)
    end if

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t *A)/(A+1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    U_n = norm2(V_n)
    V_n = V_n / U_n

    ! Sample mu and phi in CM frame
    call reac % sampleOut(mu, phi, dummy, p % E, p % pRNG)

    ! Obtain post collision speed
    V_n = rotateVector(V_n, mu, phi) * U_n

    ! Return to LAB frame
    V_n = V_n + V_cm

    ! Calculate new neutron speed and direction
    U_n = norm2(V_n)
    dir_post = V_n / U_n

    ! Update particle state and calculate mu in LAB frame
    p % E = U_n * U_n
    call p % point(dir_post)
    collDat % muL = dot_product(dir_pre, dir_post)

  end subroutine scatterFromMoving


end module neutronCEstd_class
