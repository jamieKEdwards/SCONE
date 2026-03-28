module bezierVolume_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill
  implicit none
  private

  ! Module-level counter for particles exceeding subdivision limit
  integer(shortInt), public :: bezierVolumeExceededCount = 0

  ! Maximum subdivision levels (safety limit — loop exits naturally via GJK convergence)
  integer(shortInt), parameter :: MAX_SUBDIVISIONS = 10

  ! GJK tolerance for floating-point comparisons
  real(defReal), parameter     :: GJK_EPS = 1.0E-12_defReal

  ! ---- Diagnostic: analytical sphere comparison ----
  integer(shortInt), parameter  :: DIAG_UNIT = 98
  integer(shortInt), save       :: diagTotal = 0
  integer(shortInt), save       :: diagMisclass = 0
  logical(defBool), save        :: diagFileOpen = .false.

  public :: printBezierDiagnostics

  !!
  !! 3D Bezier volume surface defined by a watertight set of bicubic Bezier patches
  !!
  !! Each patch is defined by a 4x4 grid of control points. The patches must be
  !! joined at edges in a watertight manner (shared edge control points between
  !! adjacent patches).
  !!
  !! The halfspace algorithm works as follows:
  !!   1. AABB test: if point is outside the axis-aligned bounding box -> outside
  !!   2. Global convex hull test via GJK on all control points -> outside
  !!   3. Per-patch convex hull test with uniform subdivision:
  !!      - For each patch, AABB pre-filter then GJK convex hull test on its
  !!        16 control points
  !!      - If inside any patch hull, subdivide ALL patches uniformly (watertight)
  !!      - Repeat until no patch hull contains the point
  !!   4. Ray cast against the current (watertight) corner-point mesh using
  !!      Moller-Trumbore with diverse ray directions for robustness
  !!
  !! Control points are inputted as a flat list in order:
  !!   patch1(point(1,1) point(1,2) ... point(4,4)) patch2(...) ...
  !! where each point is (x, y, z)
  !!
  !! Surface tolerance: SURF_TOL
  !!
  !! Sample dictionary input:
  !!  vol { type bezierVolume; id 1; numPatches 6; ctrlPts (x y z  x y z ... ); }
  !!
  !! Private members:
  !!   ctrlPts        -> Control points array (numPatches, 4, 4, 3)
  !!   numPatches     -> Number of Bezier patches
  !!   allPtsFlat     -> All control points as (N, 3) for global GJK test
  !!   nAllPts        -> Number of points in allPtsFlat
  !!   aabb           -> Cached axis-aligned bounding box (6)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: bezierVolume
    private
    real(defReal), dimension(:,:,:,:), allocatable :: ctrlPts
    integer(shortInt)                              :: numPatches     = 0
    ! All control points stored as (N, 3) for global GJK convex hull test
    real(defReal), dimension(:,:), allocatable     :: allPtsFlat
    integer(shortInt)                              :: nAllPts        = 0
    ! Cached AABB
    real(defReal), dimension(6)                    :: aabb = ZERO
  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
    procedure :: halfspace
    ! Internal procedures
    procedure :: inPatchAABB
    procedure :: rayCast
    procedure :: rayCastDebug
    procedure :: rayTriangle
  end type bezierVolume


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(bezierVolume), intent(in) :: self
    character(:), allocatable       :: str

    str = 'bezierVolume'

  end function myType

  !!
  !! Initialise bezierVolume from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if id < 1
  !!   fatalError if number of control points is inconsistent with numPatches
  !!
  subroutine init(self, dict)
    class(bezierVolume), intent(inout)       :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, n, m, i, j, k, l
    real(defReal), dimension(:), allocatable :: ctrlPtsList
    character(100), parameter :: Here = 'init (bezierVolume_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(ctrlPtsList, 'ctrlPts')
    call dict % get(self % numPatches, 'numPatches')

    ! Check values
    if (id < 1) then
      call fatalError(Here, 'Invalid surface id provided. ID must be >= 1')
    end if

    n = size(ctrlPtsList)

    ! Each patch has 4x4 = 16 control points, each with 3 coordinates
    if (n /= self % numPatches * 16 * 3) then
      call fatalError(Here, 'Number of control points inconsistent with numPatches. '// &
                            'Expected: '//numToChar(self % numPatches * 16 * 3)// &
                            ' Got: '//numToChar(n))
    end if

    ! Load data
    call self % setID(id)

    ! Allocate and fill control points array (numPatches, 4, 4, 3)
    allocate(self % ctrlPts(self % numPatches, 4, 4, 3))

    m = 1
    do i = 1, self % numPatches
      do j = 1, 4
        do k = 1, 4
          do l = 1, 3
            self % ctrlPts(i, j, k, l) = ctrlPtsList(m)
            m = m + 1
          end do
        end do
      end do
    end do

    ! Cache the AABB
    self % aabb = self % boundingBox()

    ! Build (N, 3) array of all control points for global GJK test
    self % nAllPts = self % numPatches * 16
    allocate(self % allPtsFlat(self % nAllPts, 3))

    m = 0
    do i = 1, self % numPatches
      do j = 1, 4
        do k = 1, 4
          m = m + 1
          self % allPtsFlat(m, :) = self % ctrlPts(i, j, k, :)
        end do
      end do
    end do

    ! Open diagnostic file
    if (.not. diagFileOpen) then
      open(unit=DIAG_UNIT, file='bezier_misclass.dat', status='replace', action='write')
      write(DIAG_UNIT, '(A)') '# x  y  z  r2  analytical_inside  bezier_inside'
      diagFileOpen = .true.
    end if

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  !! Returns bounding box enclosing all control points
  !!
  pure function boundingBox(self) result(aabb)
    class(bezierVolume), intent(in) :: self
    real(defReal), dimension(6)     :: aabb
    integer(shortInt)               :: i, j, k

    aabb(1) =  INF;  aabb(2) =  INF;  aabb(3) =  INF
    aabb(4) = -INF;  aabb(5) = -INF;  aabb(6) = -INF

    do i = 1, self % numPatches
      do j = 1, 4
        do k = 1, 4
          if (self % ctrlPts(i,j,k,1) < aabb(1)) aabb(1) = self % ctrlPts(i,j,k,1)
          if (self % ctrlPts(i,j,k,2) < aabb(2)) aabb(2) = self % ctrlPts(i,j,k,2)
          if (self % ctrlPts(i,j,k,3) < aabb(3)) aabb(3) = self % ctrlPts(i,j,k,3)
          if (self % ctrlPts(i,j,k,1) > aabb(4)) aabb(4) = self % ctrlPts(i,j,k,1)
          if (self % ctrlPts(i,j,k,2) > aabb(5)) aabb(5) = self % ctrlPts(i,j,k,2)
          if (self % ctrlPts(i,j,k,3) > aabb(6)) aabb(6) = self % ctrlPts(i,j,k,3)
        end do
      end do
    end do

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! Not available for parametric surface - returns 0
  !!
  pure function evaluate(self, r) result(c)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    c = ZERO

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! Not available yet - returns INF
  !!
  pure function distance(self, r, u) result(d)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d

    d = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! Not available - returns .false.
  !!
  pure function going(self, r, u) result(halfspace)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace

    halfspace = .false.

  end function going

  !!
  !! Check if point r is inside the AABB of a single patch (4x4 = 16 points)
  !!
  !! Args:
  !!   patch [in] -> Control points of a single patch (4, 4, 3)
  !!   r     [in] -> Point to test
  !!
  !! Result:
  !!   True if r is inside the patch AABB
  !!
  pure function inPatchAABB(self, patch, r) result(inside)
    class(bezierVolume), intent(in)                :: self
    real(defReal), dimension(4,4,3), intent(in)    :: patch
    real(defReal), dimension(3), intent(in)        :: r
    logical(defBool)                               :: inside
    real(defReal), dimension(3)                    :: minPt, maxPt
    integer(shortInt)                              :: j, k

    minPt =  INF
    maxPt = -INF

    do j = 1, 4
      do k = 1, 4
        if (patch(j,k,1) < minPt(1)) minPt(1) = patch(j,k,1)
        if (patch(j,k,2) < minPt(2)) minPt(2) = patch(j,k,2)
        if (patch(j,k,3) < minPt(3)) minPt(3) = patch(j,k,3)
        if (patch(j,k,1) > maxPt(1)) maxPt(1) = patch(j,k,1)
        if (patch(j,k,2) > maxPt(2)) maxPt(2) = patch(j,k,2)
        if (patch(j,k,3) > maxPt(3)) maxPt(3) = patch(j,k,3)
      end do
    end do

    inside = (r(1) >= minPt(1) - SURF_TOL) .and. (r(1) <= maxPt(1) + SURF_TOL) .and. &
             (r(2) >= minPt(2) - SURF_TOL) .and. (r(2) <= maxPt(2) + SURF_TOL) .and. &
             (r(3) >= minPt(3) - SURF_TOL) .and. (r(3) <= maxPt(3) + SURF_TOL)

  end function inPatchAABB

  !!
  !! Determine halfspace for a particle position using Bezier patch subdivision
  !!
  !! Algorithm:
  !!   1. AABB rejection test (cached)
  !!   2. Global GJK convex hull rejection test on all control points
  !!   3. Per-patch GJK convex hull test: for each patch, use AABB as cheap
  !!      pre-filter then GJK on the 16 control points. If inside any patch
  !!      hull, subdivide ALL patches uniformly (preserving watertightness)
  !!      and repeat.
  !!   4. When no patch hull contains the point, ray cast against the current
  !!      watertight corner-point mesh.
  !!
  !! Args:
  !!   r [in] -> Particle location
  !!   u [in] -> Particle direction
  !!
  !! Result:
  !!   True if position is in +ve halfspace (outside). False if inside.
  !!
  function halfspace(self, r, u) result(hs)
    class(bezierVolume), intent(in)                :: self
    real(defReal), dimension(3), intent(in)        :: r
    real(defReal), dimension(3), intent(in)        :: u
    logical(defBool)                               :: hs
    real(defReal), dimension(:,:,:,:), allocatable  :: currentPts
    real(defReal), dimension(16, 3)                :: patchPts
    integer(shortInt)                              :: subdivision, i, j, k, m, totalPatches
    logical(defBool)                               :: anyContaining
    ! Diagnostic locals
    real(defReal)                                  :: r2
    logical(defBool)                               :: analyticalIn, bezierIn
    real(defReal), parameter                       :: DIAG_R2 = 2.25_defReal  ! 1.5^2

    ! Step 1: AABB test
    if (r(1) < self % aabb(1) - SURF_TOL .or. r(1) > self % aabb(4) + SURF_TOL .or. &
        r(2) < self % aabb(2) - SURF_TOL .or. r(2) > self % aabb(5) + SURF_TOL .or. &
        r(3) < self % aabb(3) - SURF_TOL .or. r(3) > self % aabb(6) + SURF_TOL) then
      hs = .true.  ! outside
      goto 999
    end if

    ! Step 2: Global convex hull test via GJK
    if (.not. pointInConvexHull(self % allPtsFlat, self % nAllPts, r)) then
      hs = .true.  ! outside
      goto 999
    end if

    ! Step 3: Uniform subdivision guided by per-patch GJK hull tests
    allocate(currentPts(self % numPatches, 4, 4, 3))
    currentPts = self % ctrlPts
    totalPatches = self % numPatches

    do subdivision = 1, MAX_SUBDIVISIONS

      ! Check if point is inside any patch's convex hull
      anyContaining = .false.

      do i = 1, totalPatches
        ! Cheap AABB pre-filter
        if (.not. self % inPatchAABB(currentPts(i,:,:,:), r)) cycle

        ! GJK convex hull test on this patch's 16 control points
        m = 0
        do j = 1, 4
          do k = 1, 4
            m = m + 1
            patchPts(m, :) = currentPts(i, j, k, :)
          end do
        end do

        if (pointInConvexHull(patchPts, 16, r)) then
          anyContaining = .true.
          exit
        end if
      end do

      ! If no patch hull contains the point, go to ray cast
      if (.not. anyContaining) exit

      ! Subdivide ALL patches uniformly to maintain watertight mesh
      call subdivideAll(currentPts, totalPatches)

    end do

    ! Step 4: Ray cast against the watertight corner-point mesh
    if (subdivision > MAX_SUBDIVISIONS) then
      bezierVolumeExceededCount = bezierVolumeExceededCount + 1
    end if

    hs = .not. self % rayCast(currentPts, totalPatches, r)

    ! Debug: trace deep misclassifications
    r2 = r(1)**2 + r(2)**2 + r(3)**2
    if (r2 < 1.21_defReal .and. (.not. hs .eqv. .false.)) then
      ! Point is deep inside sphere (r < 1.1) but bezier says outside
      !$omp critical(bezierDebug)
      write(*, '(A)')           '*** DEEP MISCLASS DEBUG ***'
      write(*, '(A, 3F10.4)')   '  Point:       ', r
      write(*, '(A, F10.4)')    '  Radius:      ', sqrt(r2)
      write(*, '(A, I4)')       '  Subdiv exit: ', subdivision
      write(*, '(A, I8)')       '  Total patches:', totalPatches
      write(*, '(A, L3)')       '  hs (outside):', hs
      call self % rayCastDebug(currentPts, totalPatches, r)
      !$omp end critical(bezierDebug)
    end if

    deallocate(currentPts)

    ! ---- Diagnostic: compare with analytical sphere ----
    999 continue
    r2 = r(1)**2 + r(2)**2 + r(3)**2
    analyticalIn = (r2 < DIAG_R2)
    bezierIn     = .not. hs

    !$omp atomic
    diagTotal = diagTotal + 1

    if (analyticalIn .neqv. bezierIn) then
      !$omp critical(bezierDiag)
      diagMisclass = diagMisclass + 1
      if (diagFileOpen) then
        write(DIAG_UNIT, '(3ES16.8, ES16.8, L3, L3)') r(1), r(2), r(3), r2, analyticalIn, bezierIn
      end if
      !$omp end critical(bezierDiag)
    end if

  end function halfspace

  ! ============================================================================
  ! GJK convex hull containment test
  ! ============================================================================

  !!
  !! 3D cross product: c = a x b
  !!
  pure function cross3(a, b) result(c)
    real(defReal), dimension(3), intent(in) :: a, b
    real(defReal), dimension(3)             :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function cross3

  !!
  !! Triple cross product: (a x b) x c
  !!
  !! Useful for computing the component of c perpendicular to a, projected
  !! into the plane defined by a and b.
  !!
  pure function tripleCross(a, b, c) result(r)
    real(defReal), dimension(3), intent(in) :: a, b, c
    real(defReal), dimension(3)             :: r

    r = cross3(cross3(a, b), c)

  end function tripleCross

  !!
  !! Test if point r is inside the convex hull of a set of 3D points using GJK
  !!
  !! The Gilbert-Johnson-Keerthi algorithm tests whether the origin lies inside
  !! the convex hull of {pts(i,:) - r}, working directly on the point set without
  !! explicitly constructing the hull. Robust to degenerate configurations
  !! (coplanar, collinear, duplicate points).
  !!
  !! Args:
  !!   pts  [in] -> Array of 3D points (nPts, 3)
  !!   nPts [in] -> Number of points
  !!   r    [in] -> Query point
  !!
  !! Result:
  !!   True if r is inside the convex hull of pts
  !!
  pure function pointInConvexHull(pts, nPts, r) result(inside)
    real(defReal), dimension(:,:), intent(in) :: pts
    integer(shortInt), intent(in)             :: nPts
    real(defReal), dimension(3), intent(in)   :: r
    logical(defBool)                          :: inside
    ! GJK working variables
    real(defReal), dimension(4, 3)            :: S          ! Simplex (up to tetrahedron)
    integer(shortInt)                         :: nS         ! Current simplex size
    real(defReal), dimension(3)               :: d          ! Search direction
    real(defReal), dimension(3)               :: sup        ! Support point
    real(defReal)                             :: maxDot, dp
    integer(shortInt)                         :: i, bestIdx, iter
    integer(shortInt), parameter              :: MAX_ITER = 64

    inside = .false.

    if (nPts < 1) return

    ! Initial direction: from r toward centroid of points (in shifted space)
    d = ZERO
    do i = 1, nPts
      d = d + pts(i,:)
    end do
    d = d / real(nPts, defReal) - r

    ! If r is at the centroid, it is inside
    if (dot_product(d, d) < GJK_EPS) then
      inside = .true.
      return
    end if

    ! First support point in shifted space {pts - r}
    maxDot = -INF
    bestIdx = 1
    do i = 1, nPts
      dp = dot_product(pts(i,:) - r, d)
      if (dp > maxDot) then
        maxDot = dp
        bestIdx = i
      end if
    end do
    sup = pts(bestIdx,:) - r

    ! If support doesn't reach origin along d, origin is outside
    if (dot_product(sup, d) < ZERO) return

    ! Initialise simplex with first support point
    nS = 1
    S(1,:) = sup
    d = -sup  ! Direction from support toward origin

    ! Main GJK loop
    do iter = 1, MAX_ITER
      ! Find support point in direction d
      maxDot = -INF
      bestIdx = 1
      do i = 1, nPts
        dp = dot_product(pts(i,:) - r, d)
        if (dp > maxDot) then
          maxDot = dp
          bestIdx = i
        end if
      end do
      sup = pts(bestIdx,:) - r

      ! If new support doesn't pass origin along d, origin is outside hull
      if (dot_product(sup, d) < ZERO) return

      ! Add to simplex
      nS = nS + 1
      S(nS,:) = sup

      ! Process simplex — updates S, nS, d, and may set inside = .true.
      call gjkDoSimplex(S, nS, d, inside)
      if (inside) return
    end do

    ! Failed to converge — conservatively report inside (triggers subdivision)
    inside = .true.

  end function pointInConvexHull

  !!
  !! Process the GJK simplex: determine new search direction or detect containment
  !!
  !! The newest point A is always at position nS. Handles line (nS=2),
  !! triangle (nS=3), and tetrahedron (nS=4) cases. May reduce the simplex
  !! by discarding vertices that are not closest to the origin.
  !!
  !! Args:
  !!   S      [inout] -> Simplex vertices, may be rearranged/reduced
  !!   nS     [inout] -> Simplex size, may decrease
  !!   d      [inout] -> New search direction toward origin
  !!   inside [out]   -> True if origin is inside the simplex (tetrahedron case)
  !!
  pure subroutine gjkDoSimplex(S, nS, d, inside)
    real(defReal), dimension(4,3), intent(inout) :: S
    integer(shortInt), intent(inout)             :: nS
    real(defReal), dimension(3), intent(inout)   :: d
    logical(defBool), intent(out)                :: inside
    real(defReal), dimension(3) :: a, b, c, dd
    real(defReal), dimension(3) :: ab, ac, ad, ao
    real(defReal), dimension(3) :: abc, acd, adb
    real(defReal), dimension(3) :: abPerp, acPerp, tcross
    logical(defBool)            :: checkAB

    inside = .false.

    ! --- Tetrahedron case (nS = 4) ---
    ! Check which face the origin is outside; reduce to that triangle.
    ! If inside all faces -> origin is inside the tetrahedron.
    if (nS == 4) then
      a  = S(4,:)
      b  = S(3,:)
      c  = S(2,:)
      dd = S(1,:)
      ab = b - a;  ac = c - a;  ad = dd - a;  ao = -a

      ! Face normals, oriented outward from the tetrahedron
      abc = cross3(ab, ac)
      if (dot_product(abc, ad) > ZERO) abc = -abc   ! away from D
      acd = cross3(ac, ad)
      if (dot_product(acd, ab) > ZERO) acd = -acd   ! away from B
      adb = cross3(ad, ab)
      if (dot_product(adb, ac) > ZERO) adb = -adb   ! away from C

      if (dot_product(abc, ao) > ZERO) then
        ! Outside face ABC — reduce to triangle, preserve outward winding
        if (dot_product(cross3(ab, ac), ad) < ZERO) then
          S(1,:) = c;  S(2,:) = b;  S(3,:) = a
        else
          S(1,:) = b;  S(2,:) = c;  S(3,:) = a
        end if
        nS = 3

      else if (dot_product(acd, ao) > ZERO) then
        ! Outside face ACD
        if (dot_product(cross3(ac, ad), ab) < ZERO) then
          S(1,:) = dd;  S(2,:) = c;  S(3,:) = a
        else
          S(1,:) = c;  S(2,:) = dd;  S(3,:) = a
        end if
        nS = 3

      else if (dot_product(adb, ao) > ZERO) then
        ! Outside face ADB
        if (dot_product(cross3(ad, ab), ac) < ZERO) then
          S(1,:) = b;  S(2,:) = dd;  S(3,:) = a
        else
          S(1,:) = dd;  S(2,:) = b;  S(3,:) = a
        end if
        nS = 3

      else
        ! Inside all faces — origin is enclosed
        inside = .true.
        return
      end if
    end if

    ! --- Triangle case (nS = 3) ---
    ! Determine whether origin is closest to edge AC, edge AB, or the face.
    if (nS == 3) then
      a  = S(3,:)
      b  = S(2,:)
      c  = S(1,:)
      ab = b - a;  ac = c - a;  ao = -a

      abc    = cross3(ab, ac)         ! face normal
      acPerp = cross3(abc, ac)        ! perpendicular to AC, pointing away from B
      abPerp = cross3(ab, abc)        ! perpendicular to AB, pointing away from C

      checkAB = .false.

      if (dot_product(acPerp, ao) > ZERO) then
        ! Outside edge AC
        if (dot_product(ac, ao) > ZERO) then
          ! Closest to edge AC
          S(1,:) = c;  S(2,:) = a
          nS = 2
          d = tripleCross(ac, ao, ac)
          return
        else
          checkAB = .true.
        end if
      end if

      if (checkAB .or. dot_product(abPerp, ao) > ZERO) then
        ! Outside edge AB or vertex A region
        if (dot_product(ab, ao) > ZERO) then
          S(1,:) = b;  S(2,:) = a
          nS = 2
          d = tripleCross(ab, ao, ab)
        else
          S(1,:) = a
          nS = 1
          d = ao
        end if
        return
      end if

      ! Inside the triangle — origin is above or below the face
      if (dot_product(abc, ao) > ZERO) then
        d = abc
      else
        ! Flip winding so normal points toward origin
        S(1,:) = b;  S(2,:) = c;  S(3,:) = a
        d = -abc
      end if
      return
    end if

    ! --- Line case (nS = 2) ---
    if (nS == 2) then
      a  = S(2,:)
      b  = S(1,:)
      ab = b - a;  ao = -a

      if (dot_product(ab, ao) > ZERO) then
        ! Origin is alongside the edge — search perpendicular to AB toward origin
        tcross = tripleCross(ab, ao, ab)
        if (dot_product(tcross, tcross) < GJK_EPS) then
          ! AB and AO are parallel — origin is on the line segment
          ! Pick any perpendicular direction
          if (abs(ab(1)) < abs(ab(2))) then
            d(1) = ZERO;  d(2) = -ab(3);  d(3) = ab(2)
          else
            d(1) = ab(3);  d(2) = ZERO;  d(3) = -ab(1)
          end if
        else
          d = tcross
        end if
      else
        ! Origin is behind A — reduce to just A
        S(1,:) = a
        nS = 1
        d = ao
      end if
      return
    end if

  end subroutine gjkDoSimplex

  ! ============================================================================
  ! Bezier patch subdivision
  ! ============================================================================

  !!
  !! Subdivide all patches into 4 sub-patches each using De Casteljau algorithm
  !!
  !! Each bicubic patch is subdivided in both u and v at t=0.5, producing 4 sub-patches.
  !! The output array is reallocated to hold 4x as many patches.
  !!
  !! Args:
  !!   pts      [inout] -> Control points, reallocated to (4*nPatches, 4, 4, 3)
  !!   nPatches [inout] -> Updated to 4*nPatches
  !!
  subroutine subdivideAll(pts, nPatches)
    real(defReal), dimension(:,:,:,:), allocatable, intent(inout) :: pts
    integer(shortInt), intent(inout)                          :: nPatches
    real(defReal), dimension(:,:,:,:), allocatable            :: newPts
    real(defReal), dimension(4,4,3)                           :: patch
    real(defReal), dimension(4,4,3)                           :: Q00, Q01, Q10, Q11
    integer(shortInt)                                         :: i, newN

    newN = nPatches * 4
    allocate(newPts(newN, 4, 4, 3))

    do i = 1, nPatches
      patch = pts(i,:,:,:)
      call subdividePatch(patch, Q00, Q01, Q10, Q11)
      newPts((i-1)*4 + 1, :,:,:) = Q00
      newPts((i-1)*4 + 2, :,:,:) = Q01
      newPts((i-1)*4 + 3, :,:,:) = Q10
      newPts((i-1)*4 + 4, :,:,:) = Q11
    end do

    call move_alloc(newPts, pts)
    nPatches = newN

  end subroutine subdivideAll

  !!
  !! Subdivide a single bicubic Bezier patch into 4 sub-patches at (u,v) = (0.5, 0.5)
  !!
  !! Uses De Casteljau algorithm, first in u then in v direction.
  !!
  !! Args:
  !!   patch [in]  -> Original 4x4 control points
  !!   Q00   [out] -> Bottom-left sub-patch
  !!   Q01   [out] -> Bottom-right sub-patch
  !!   Q10   [out] -> Top-left sub-patch
  !!   Q11   [out] -> Top-right sub-patch
  !!
  subroutine subdividePatch(patch, Q00, Q01, Q10, Q11)
    real(defReal), dimension(4,4,3), intent(in)  :: patch
    real(defReal), dimension(4,4,3), intent(out) :: Q00, Q01, Q10, Q11
    real(defReal), dimension(4,4,3)              :: left, right
    real(defReal), dimension(4,3)                :: rowL, rowR
    real(defReal), dimension(4,3)                :: colL, colR
    integer(shortInt)                            :: i, j

    ! First subdivide each row in u direction at t=0.5
    do i = 1, 4
      call subdivideRow(patch(i,:,:), rowL, rowR)
      left(i,:,:)  = rowL
      right(i,:,:) = rowR
    end do

    ! Then subdivide each column of left and right halves in v direction at t=0.5
    do j = 1, 4
      call subdivideRow(left(:,j,:), colL, colR)
      Q00(:,j,:) = colL
      Q10(:,j,:) = colR
    end do

    do j = 1, 4
      call subdivideRow(right(:,j,:), colL, colR)
      Q01(:,j,:) = colL
      Q11(:,j,:) = colR
    end do

  end subroutine subdividePatch

  !!
  !! Subdivide a cubic Bezier row of 4 control points at t=0.5 using De Casteljau
  !!
  !! Args:
  !!   row  [in]  -> 4 control points (4, 3)
  !!   left [out] -> Left sub-curve control points
  !!   right[out] -> Right sub-curve control points
  !!
  subroutine subdivideRow(row, left, right)
    real(defReal), dimension(4,3), intent(in)  :: row
    real(defReal), dimension(4,3), intent(out) :: left, right
    real(defReal), dimension(3)                :: p01, p12, p23, p012, p123, p0123

    p01   = HALF * (row(1,:) + row(2,:))
    p12   = HALF * (row(2,:) + row(3,:))
    p23   = HALF * (row(3,:) + row(4,:))
    p012  = HALF * (p01  + p12)
    p123  = HALF * (p12  + p23)
    p0123 = HALF * (p012 + p123)

    left(1,:)  = row(1,:)
    left(2,:)  = p01
    left(3,:)  = p012
    left(4,:)  = p0123

    right(1,:) = p0123
    right(2,:) = p123
    right(3,:) = p23
    right(4,:) = row(4,:)

  end subroutine subdivideRow

  ! ============================================================================
  ! Ray casting
  ! ============================================================================

  !!
  !! Ray cast against corner-point mesh to determine inside/outside
  !!
  !! Fires a ray from r and counts intersections with the triangulated corner
  !! mesh (2 triangles per patch quad). Odd intersections -> inside, Even -> outside.
  !! If degenerate intersection detected, perturbs ray direction and retries.
  !!
  !! Args:
  !!   pts      [in] -> Current control points
  !!   nPatches [in] -> Number of patches
  !!   r        [in] -> Point to test
  !!
  !! Result:
  !!   True if inside
  !!
  function rayCast(self, pts, nPatches, r) result(inside)
    class(bezierVolume), intent(in)               :: self
    real(defReal), dimension(:,:,:,:), intent(in) :: pts
    integer(shortInt), intent(in)                 :: nPatches
    real(defReal), dimension(3), intent(in)       :: r
    logical(defBool)                              :: inside
    real(defReal), dimension(3)                   :: rayDir
    real(defReal), dimension(3)                   :: v0, v1, v2, v3
    integer(shortInt)                             :: count, i
    logical(defBool)                              :: hit, nearZero, anyNearZero
    real(defReal), dimension(3)                   :: rLocal
    real(defReal), parameter                      :: NUDGE = 1.0E-5_defReal

    rayDir = (/ ONE, ONE / 3.0_defReal, ONE / 7.0_defReal /)
    rayDir = rayDir / norm2(rayDir)

    ! First pass with original point
    count = 0
    anyNearZero = .false.

    do i = 1, nPatches
      v0 = pts(i, 1, 1, :)
      v1 = pts(i, 1, 4, :)
      v2 = pts(i, 4, 4, :)
      v3 = pts(i, 4, 1, :)

      call self % rayTriangle(r, rayDir, v0, v1, v2, hit, nearZero)
      if (hit) count = count + 1
      if (nearZero) anyNearZero = .true.

      call self % rayTriangle(r, rayDir, v0, v2, v3, hit, nearZero)
      if (hit) count = count + 1
      if (nearZero) anyNearZero = .true.
    end do

    ! If degenerate (count=0 and near-zero detected), nudge and recast once
    if (count == 0 .and. anyNearZero) then
      rLocal = r + (/ NUDGE, NUDGE, NUDGE /)
      count = 0
      anyNearZero = .false.
      do i = 1, nPatches
        v0 = pts(i, 1, 1, :)
        v1 = pts(i, 1, 4, :)
        v2 = pts(i, 4, 4, :)
        v3 = pts(i, 4, 1, :)

        call self % rayTriangle(rLocal, rayDir, v0, v1, v2, hit, nearZero)
        if (hit) count = count + 1
        if (nearZero) anyNearZero = .true.

        call self % rayTriangle(rLocal, rayDir, v0, v2, v3, hit, nearZero)
        if (hit) count = count + 1
        if (nearZero) anyNearZero = .true.
      end do
    end if

    inside = mod(count, 2) == 1

  end function rayCast

  !!
  !! Debug version of rayCast — prints per-patch hit details
  !!
  subroutine rayCastDebug(self, pts, nPatches, r)
    class(bezierVolume), intent(in)               :: self
    real(defReal), dimension(:,:,:,:), intent(in) :: pts
    integer(shortInt), intent(in)                 :: nPatches
    real(defReal), dimension(3), intent(in)       :: r
    real(defReal), dimension(3)                   :: rayDir
    real(defReal), dimension(3)                   :: v0, v1, v2, v3
    integer(shortInt)                             :: count, i
    logical(defBool)                              :: hit1, hit2, nz

    rayDir = (/ ONE, ONE / 3.0_defReal, ONE / 7.0_defReal /)
    rayDir = rayDir / norm2(rayDir)

    count = 0
    write(*, '(A, 3F10.6)') '  Ray direction: ', rayDir

    do i = 1, nPatches
      v0 = pts(i, 1, 1, :)
      v1 = pts(i, 1, 4, :)
      v2 = pts(i, 4, 4, :)
      v3 = pts(i, 4, 1, :)

      call self % rayTriangle(r, rayDir, v0, v1, v2, hit1, nz)
      call self % rayTriangle(r, rayDir, v0, v2, v3, hit2, nz)

      if (hit1 .or. hit2) then
        count = count + merge(1, 0, hit1) + merge(1, 0, hit2)
        write(*, '(A, I6, A, L2, A, L2)') '  Patch ', i, ': T1=', hit1, ' T2=', hit2
        write(*, '(A, 3F10.4)') '    v0=', v0
        write(*, '(A, 3F10.4)') '    v1=', v1
        write(*, '(A, 3F10.4)') '    v2=', v2
        write(*, '(A, 3F10.4)') '    v3=', v3
      end if
    end do

    write(*, '(A, I4, A, L3)') '  Total hits: ', count, '  inside=', mod(count, 2) == 1

  end subroutine rayCastDebug

  !!
  !! Watertight ray-triangle intersection (Woop, Benthin & Wald 2013)
  !!
  !! Transforms triangle vertices into ray-aligned coordinates and computes
  !! edge functions that depend only on each edge's two vertices. Shared edges
  !! between adjacent triangles produce bitwise-identical edge function values,
  !! guaranteeing exactly one hit per shared edge regardless of ray direction.
  !!
  !! Args:
  !!   r          [in]  -> Ray origin
  !!   rayDir     [in]  -> Ray direction (normalised)
  !!   v0,v1,v2   [in]  -> Triangle vertices
  !!   hit        [out] -> True if ray intersects triangle at t > 0
  !!   nearZero   [out] -> True if ray parameter t/det is near zero (point on triangle plane)
  !!
  subroutine rayTriangle(self, r, rayDir, v0, v1, v2, hit, nearZero)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, rayDir, v0, v1, v2
    logical(defBool), intent(out)           :: hit, nearZero
    real(defReal), dimension(3)             :: A, B, C
    real(defReal)                           :: Ax, Ay, Bx, By, Cx, Cy
    real(defReal)                           :: Az, Bz, Cz
    real(defReal)                           :: e0, e1, e2, det, t
    real(defReal)                           :: absDir1, absDir2, absDir3
    integer(shortInt)                       :: kz, kx, ky
    real(defReal)                           :: Sx, Sy, Sz
    real(defReal), parameter                :: EPS = 1.0E-10_defReal
    real(defReal), parameter                :: NEAR_ZERO_TOL = 1.0E-4_defReal

    hit = .false.
    nearZero = .false.

    ! Step 1: Find dominant axis of ray direction for coordinate permutation
    absDir1 = abs(rayDir(1))
    absDir2 = abs(rayDir(2))
    absDir3 = abs(rayDir(3))

    if (absDir1 > absDir2 .and. absDir1 > absDir3) then
      kz = 1; kx = 2; ky = 3
    else if (absDir2 > absDir3) then
      kz = 2; kx = 3; ky = 1
    else
      kz = 3; kx = 1; ky = 2
    end if

    ! Step 2: Shear constants to align ray with +z axis
    Sx = rayDir(kx) / rayDir(kz)
    Sy = rayDir(ky) / rayDir(kz)
    Sz = ONE / rayDir(kz)

    ! Step 3: Translate vertices relative to ray origin
    A = v0 - r
    B = v1 - r
    C = v2 - r

    ! Step 4: Shear and permute to ray-aligned coordinates
    Ax = A(kx) - Sx * A(kz)
    Ay = A(ky) - Sy * A(kz)
    Bx = B(kx) - Sx * B(kz)
    By = B(ky) - Sy * B(kz)
    Cx = C(kx) - Sx * C(kz)
    Cy = C(ky) - Sy * C(kz)

    ! Step 5: Edge functions — each depends only on its two vertices
    ! This is the key to watertightness: shared edges give identical results
    e0 = Bx * Cy - By * Cx   ! edge v1-v2
    e1 = Cx * Ay - Cy * Ax   ! edge v2-v0
    e2 = Ax * By - Ay * Bx   ! edge v0-v1

    ! Step 6: Check if point is inside triangle
    ! All edge functions must have the same sign (or be zero)
    ! Use a small tolerance so that near-zero edge values (point on/near edge)
    ! are treated as zero — prevents both adjacent triangles rejecting the point
    if (e0 < -EPS .or. e1 < -EPS .or. e2 < -EPS) then
      if (e0 > EPS .or. e1 > EPS .or. e2 > EPS) return
    end if

    ! Determinant
    det = e0 + e1 + e2
    if (abs(det) < EPS) then
      nearZero = .true.
      return
    end if

    ! Step 7: Compute t (ray parameter) — only z-shear needed now
    Az = Sz * A(kz)
    Bz = Sz * B(kz)
    Cz = Sz * C(kz)
    t = e0 * Az + e1 * Bz + e2 * Cz

    ! Flag near-zero t: point lies on or very near the triangle plane
    ! This is checked BEFORE the forward-direction filter so we detect
    ! rejected-but-on-plane cases that cause the 0-hit degeneracy
    if (abs(t) < NEAR_ZERO_TOL * abs(det)) nearZero = .true.

    ! Check forward intersection (accounting for sign of det)
    if (det > ZERO) then
      if (t < EPS * det) return
    else
      if (t > EPS * det) return
    end if

    hit = .true.

  end subroutine rayTriangle

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(bezierVolume), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    if (allocated(self % ctrlPts))    deallocate(self % ctrlPts)
    if (allocated(self % allPtsFlat)) deallocate(self % allPtsFlat)
    self % numPatches = 0
    self % nAllPts    = 0
    self % aabb       = ZERO

  end subroutine kill

  !!
  !! Print diagnostic summary to stdout and close file
  !! Call this at the end of the run to get the misclassification report
  !!
  subroutine printBezierDiagnostics()
    real(defReal) :: pct

    if (diagTotal == 0) return

    if (diagTotal > 0) then
      pct = 100.0_defReal * real(diagMisclass, defReal) / real(diagTotal, defReal)
    else
      pct = ZERO
    end if

    write(*, '(A)')          ''
    write(*, '(A)')          '====== bezierVolume Diagnostic Report ======'
    write(*, '(A, I12)')     'Total halfspace calls:  ', diagTotal
    write(*, '(A, I12)')     'Misclassified points:   ', diagMisclass
    write(*, '(A, F10.4, A)') 'Misclassification rate: ', pct, '%'
    write(*, '(A)')          'Details written to: bezier_misclass.dat'
    write(*, '(A)')          '============================================'

    ! Write summary to file and close
    if (diagFileOpen) then
      write(DIAG_UNIT, '(A)')          '# ---- Summary ----'
      write(DIAG_UNIT, '(A, I12)')     '# Total calls:    ', diagTotal
      write(DIAG_UNIT, '(A, I12)')     '# Misclassified:  ', diagMisclass
      write(DIAG_UNIT, '(A, F10.4, A)') '# Rate:           ', pct, '%'
      close(DIAG_UNIT)
      diagFileOpen = .false.
    end if
  end subroutine printBezierDiagnostics

end module bezierVolume_class
