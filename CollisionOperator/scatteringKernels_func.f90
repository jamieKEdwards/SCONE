module scatteringKernels_func

  use numPrecision
  use genericProcedures,          only : rotateVector
  use RNG_class,                  only : RNG
  use particle_class,             only : particle
  use aceNeutronNuclide_class,    only : aceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,      only : aceNeutronNuclide

  implicit none
  private

  public  :: asymptoticScatter
  public  :: targetVelocity_constXS
  public  :: targetVelocity_DBRCXS
  public  :: asymptoticInelasticScatter

  private :: sample_x2expx2
  private :: sample_x3expx2

  !type(aceNeutronNuclide), pointer, public :: aceNuc => null()


  !type(aceNeutronNuclide),dimension(:),pointer :: nuclides  => null()

contains


  !!
  !! Subroutine to perform ELASTIC SCATTERING from a stationary target nucleus.
  !! Changes direction and energy in LAB system
  !! Takes mu in CM system, returns it in LAB system
  !!
  !! Dir needs to be normalised. Prodedure will produce wrong results without error message
  !! if it is not.
  !!
  !! Based on MCNP manual chapter 2.
  !!
  subroutine asymptoticScatter(E,mu,A)
    real(defReal), intent(inout)               :: E    !! Pre-collision energy in LAB
    real(defReal), intent(inout)               :: mu   !! Cosine of delection angle
    real(defReal), intent(in)                  :: A    !! Target mass [neutrons]
    real(defReal)                              :: E_in
    real(defReal)                              :: inv_Ap1

    ! Store initial energy and precalculate 1/(A+1)
    E_in = E
    inv_Ap1 = 1.0/ (A + 1.0)

    ! Find post-collision energy
    E  = (1.0 + A*A + 2 *A*mu) *E_in * inv_Ap1 * inv_Ap1

    ! Find deflection angle in LAB
    mu = (A*mu + 1)*sqrt(E_in/E)* inv_Ap1

    ! Correct possible nuclear data shortcomings
    if (mu > ONE) mu = ONE

  end subroutine asymptoticScatter

  !!
  !! Subroutine to perform INELASTIC SCATTERING from a stationary target nucleus.
  !! Changes direction and energy in LAB system
  !! Takes mu in CM system, and E_out in CM system.
  !! Returns mu in LAB system.
  !!
  !! Dir needs to be normalised. Prodedure will produce wrong results without error message
  !! if it is not.
  !!
  !! Based on MCNP manual chapter 2.
  !!
  subroutine asymptoticInelasticScatter(E,mu,E_out,A)
    real(defReal), intent(inout) :: E     !! Pre-collision energy in Lab
    real(defReal), intent(inout) :: mu    !! Cosine of delection angle
    real(defReal), intent(in)    :: E_out !! Post collision energy in CM
    real(defReal), intent(in)    :: A     !! Target mass [neutrons]
    real(defReal)                :: E_in
    real(defReal)                :: inv_Ap1

    ! Store initial energy and precalculate 1/(A+1)
    E_in = E
    inv_Ap1 = 1.0/ (A + 1.0)

    ! Find post-collision energy
    E = E_out + (E_in +TWO*mu*(A+ONE)*sqrt(E_in*E_out)) * inv_Ap1 * inv_Ap1

    ! Find deflection angle in LAB
    mu = mu * sqrt(E_out/E) + sqrt(E_in/E)* inv_Ap1

    ! Correct possible nuclear data shortcomings
    if (mu > ONE) mu = ONE

  end subroutine asymptoticInelasticScatter


  !!
  !! Function that returns a sample of target velocity using constant XS approximation
  !! V_t is a vector. The velocity is scaled by a factor sqrt(Mn/2) where Mn is mass of a neutron
  !! so that V_t*V_t=E_t with E_t beeing kinetic energy of a NEUTRON traveling with TARGET VELOCITY
  !! (note that it is not a kinetic energy of the target).
  !!
  !!
  function targetVelocity_constXS(E,dir,A,kT,rand) result (V_t)
    real(defReal), intent(in)               :: E
    real(defReal), dimension(3), intent(in) :: dir
    real(defReal), intent(in)               :: A
    real(defReal), intent(in)               :: kT
    class(RNG), intent(inout)               :: rand
    real(defReal),dimension(3)              :: V_t
    real(defReal)                           :: alpha, mu, phi, P_acc
    real(defReal)                           :: X, Y
    real(defReal)                           :: r1, r2, r3, r4

    ! Calculate neutron Y = beta *V_n
    ! beta = sqrt(A*Mn/2kT). Note velocity scaling by sqrt(Mn/2).
    Y = sqrt(A*E/kT)

    ! Calculate treshhold factor alpha
    alpha = 2.0/(Y*sqrt(PI)+2.0)


    rejectionLoop: do

      ! Obtain random numbers
      r1 = rand % get()
      r2 = rand % get()
      r3 = rand % get()

      ! Sample X = beta * V_t
      if ( r1 > alpha ) then
        X = sample_x2expx2(rand)

      else
        X = sample_x3expx2(rand)

      end if

      ! Sample polar angle of target velocity wrt. neutron direction
      mu = 2.0 * r2 - 1.0;

      ! Calculate Acceptance Propability
      P_acc = sqrt(Y*Y + X*X - 2.0*X*Y*mu) / (Y+X)

      ! Accept or reject mu
      if (P_acc > r3) exit rejectionLoop

    end do rejectionLoop

    ! Calculate azimithal angle for traget and obtain target direction
    r4 = rand % get()
    phi = 2.0 * PI * r4

    V_t = rotateVector(dir, mu, phi)

    ! Scale target direction by magnitude of velocity
    V_t = V_t * (X * sqrt(kT/A))

  end function targetVelocity_constXS



  !!
  !! Function that returns a sample of target velocity with DBRC
  !! V_t is a vector. The velocity is scaled by a factor sqrt(Mn/2) where Mn is mass of a neutron
  !! so that V_t*V_t=E_t with E_t being kinetic energy of a NEUTRON traveling with TARGET VELOCITY
  !! (note that it is not a kinetic energy of the target).
  !!
  !!
  function targetVelocity_DBRCXS(aceNuc, E, dir, A, kT, rand, TmajXS) result (V_t)
    class(aceNeutronNuclide), intent(in)    :: aceNuc
    real(defReal), intent(in)               :: E
    real(defReal), dimension(3), intent(in) :: dir
    real(defReal), intent(in)               :: A
    real(defReal), intent(in)               :: kT
    real(defReal), intent(in)               :: TmajXS
    !integer(shortInt), intent(in)           :: nucIdx
    class(RNG), intent(inout)               :: rand
    !class(aceNeutronNuclide), pointer       :: ptr
    integer(shortInt)                       :: idx
    real(defReal),dimension(3)              :: V_t
    real(defReal)                           :: alpha, mu, phi, P_acc, DBRC_acc
    real(defReal)                           :: X, Y
    real(defReal)                           :: r1, r2, r3, r4, r5
    real(defreal)                           :: rel_v, rel_E, xs_rel_v, f

    !print *, "DBRC Target velocity function being used"
    ! Pointer to aceNeutronNuclide
    !associate(nuc => aceNeutronNuclide(nucIdx))

    ! Calculate neutron Y = beta *V_n
    ! beta = sqrt(A*Mn/2kT). Note velocity scaling by sqrt(Mn/2).
    Y = sqrt(A * E / kT)

    ! Calculate treshhold factor alpha
    ! In MCNP, alpha is p1
    alpha = 2.0 / (Y * sqrt(PI) + 2.0)

    ! Call through temperature majorant
    !call xsData % updateTempMajorantXS(E, kT, A)

    rejectionLoop: do
      !print*, "Rejected"
      ! Obtain random numbers
      r1 = rand % get()
      r2 = rand % get()
      r3 = rand % get()
      r4 = rand % get()

      ! Sample X = beta * V_t
      if ( r1 > alpha ) then
        X = sample_x2expx2(rand)

      else
        X = sample_x3expx2(rand)

      end if

      ! Sample polar angle of target velocity wrt. neutron direction
      mu = 2.0 * r2 - 1.0;

      ! Calculate relative velocity between neutron and target
      rel_v = sqrt(Y * Y + X * X - 2.0 * X * Y * mu)

      ! Calculate Acceptance Propability
      P_acc = rel_v / (Y + X)

      ! Accept or reject mu
      if (P_acc < r3) cycle

      !print*, rel_v
      ! Relative energy = relative velocity **2 due to sqrt(Mn/2) scaling factor
      rel_E = (rel_v**2 * kT / A)


      !print*, rel_E
      ! Find scattering xs of target at relative velocity (rel_v), done through calling subroutines from neutronNuclide
      call aceNuc % search(idx, f, rel_E)
      !call aceNuc % scatterXS(xs, idx, f)
      xs_rel_v = aceNuc % scatterXS(idx, f)
      !print *, "rel micro scatter xs", xs_rel_v

      ! Introduce DBRC acceptance condition
      ! first term is ratio of cross sections, second term is the C' term from Becker 2009
      DBRC_acc = (xs_rel_v / TmajXS) * ((1 + (Y * sqrt(PI)) / 2) / (Y * sqrt(PI)))
      !print *, "DBRC_acc", DBRC_acc

      ! accept or reject with DBRC
      if (DBRC_acc > r4) then
        !print *, "ACCEPTED"
        exit rejectionLoop
      else
        !print *, "rejected"
      end if


    end do rejectionLoop

    !End associate
    !print*, "Accepted"
    ! Calculate azimithal angle for traget and obtain target direction
    r5 = rand % get()
    phi = 2.0 * PI * r5

    V_t = rotateVector(dir, mu, phi)

    ! Scale target direction by magnitude of velocity
    V_t = V_t * (X * sqrt(kT / A))

  end function targetVelocity_DBRCXS








  !!
  !! Helper function to sample x^2 * exp( - x^2) probability distribution
  !! Uses random numbers from provided RNG
  !! Changes variables to y = x^2 which transforms PDF to Gamma(3/2,1) distribution
  !! Then uses method based on Johnk's theorem and sum of Gamma distributed random variables
  !!
  function sample_x2expx2(rand) result(sample)
    class(RNG), intent(inout) :: rand
    real(defReal)             :: sample
    real(defReal)             :: r1, r2, r3
    real(defReal)             :: beta, gamma05, cosine

    ! Obtain all random numbers
    r1 = rand % get()
    r2 = rand % get()
    r3 = rand % get()

    ! Define B(a,b) as a RANDOM VARIABLE governed by beta(a,b) distribution
    ! Similarly define G(a,b) as a RANDOM VARIABLE governed by gamma(a,b) distribution
    ! Note that B(a,b)*G(a,b) is a product of RANDOM VARIABLES not of Probability Density Functions

    ! Obtain sample of B(0.5,0.5) distribution based on Johnk's Theorem. Sample of cosine
    ! instead of using rejection scheme
    cosine = cos(0.5*PI*r1)
    beta = cosine * cosine

    ! Obtain sample of G(0.5,1) using the fact that G(0.5,1) = B(0.5,0.5) * G(1,1)
    ! G(1,1) is just exponential distribution
    gamma05 = -log(r2) * beta

    ! Obtain sample of G(3/2,1) using the facte that G(3/2,1) = G(0.5,1) + G(1,1)
    sample = -log(r3) + gamma05

    ! Change variables back to x from y=x^2
    sample = sqrt(sample)

  end function sample_x2expx2


  !!
  !! Helper function to sample x^3 * exp( - x^2) probability distribution
  !! Uses random numbers from provided RNG
  !! Changes variables to y = x^2 which transforms PDF to Gamma(2,1) distribution
  !! Sampling Gamma(2,1) is trivial using sum of Gamma distributed random variables
  !!
  function sample_x3expx2(rand) result(sample)
    class(RNG), intent(inout) :: rand
    real(defReal)             :: sample
    real(defReal)             :: r1, r2

    ! Obtain random numbers
    r1 = rand % get()
    r2 = rand % get()

    ! Sample Gamma(2,1) by summing two samples of Gamma(1,1) [exponential distribution]
    sample = -log(r1) - log(r2)

    ! Change variables back to x
    sample = sqrt(sample)

  end function sample_x3expx2


end module scatteringKernels_func
