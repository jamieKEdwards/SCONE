module bezierVolume_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL, NUDGE
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill
  use bezierPatch_func,   only : subdividePatch
  implicit none
  private

  ! Safety limit on subdivision depth per halfspace call
  integer(shortInt), parameter :: MAX_SUBDIVISIONS = 10

  ! Max sub-patch corners collected along one edge for T-junction stitching.
  ! One side subdivided to the full MAX_SUBDIVISIONS depth holds at most
  ! 2**MAX_SUBDIVISIONS + 1 distinct points; mergeEdgeVerts can merge a second,
  ! independently-subdivided side into the same buffer, so this covers double
  ! that worst case. insertUnique fatalErrors if this is ever still exceeded.
  integer(shortInt), parameter :: MAX_TJ_BUF = 2 * (2**MAX_SUBDIVISIONS) + 2

  ! Polygon buffer capacity for fanTriangulate's fan triangulation of each
  ! triangle face -- two edge buffers concatenated, so sized off MAX_TJ_BUF.
  integer(shortInt), parameter :: MAX_POLY = 2 * MAX_TJ_BUF + 1

  !!
  !! The current (possibly non-uniformly subdivided) sub-patch list produced
  !! by subdivide and threaded, unmodified, through the rest of
  !! the halfspace mesh-building pipeline.
  !!
  !! Public Members:
  !!   ctrlPts -> sub-patch control points, ctrlPts(1:n,4,4,3)
  !!   orig    -> which original patch each sub-patch descends from, orig(1:n)
  !!   uvRange -> each sub-patch's UV range within its original patch, uvRange(1:n,4)
  !!   n       -> number of valid sub-patches in ctrlPts/orig/uvRange
  !!
  type :: subPatches
    real(defReal), dimension(:,:,:,:), allocatable :: ctrlPts
    integer(shortInt), dimension(:), allocatable    :: orig
    real(defReal), dimension(:,:), allocatable      :: uvRange
    integer(shortInt) :: n = 0
  end type subPatches

  !!
  !! 3D Bezier volume surface defined by a watertight set of bicubic Bezier patches.
  !!
  !! The closed solid enclosed by this set of patches is sometimes called a
  !! "Bezierhedron" -- by analogy with a polyhedron, but with curved Bezier
  !! faces instead of flat polygonal ones.
  !!
  !! Each patch is a 4x4 grid of control points. Adjacent patches must share their
  !! boundary control points exactly (watertight input mesh).
  !!
  !! Halfspace algorithm:
  !!   1. AABB rejection: outside bounding box -> outside
  !!   2. Global GJK rejection: outside convex hull of all control points -> outside
  !!   3. Selective subdivision: only subdivide patches whose GJK hull contains the
  !!      query point. Track which original patches were subdivided. If no patch's
  !!      hull contains the point (typically floating-point precision on an
  !!      already-tight, deeply-subdivided hull, near a shared seam), fall back to
  !!      AABB-only selection, scoped to patches already hull-selected this call or
  !!      adjacent to one (everHullSelected) so unrelated patches are not pulled in
  !!      by loose bounding-box coincidence.
  !!   4. Build watertight triangle mesh via polygon-fan triangulation:
  !!      - Subdivided patches: fan each sub-patch, collecting T-junction vertices
  !!        both within the same original patch (collectIntraPatchEdgeVerts) and,
  !!        along the original patch's true outer boundary, from an adjacent patch
  !!        that was ALSO independently subdivided (mergeEdgeVerts) -- two adjacent
  !!        patches subdivided under the same query point are not otherwise
  !!        synchronised, so without this their independent splits can leave a gap.
  !!      - Unsubdivided patches: fan from C00, collecting T-junction vertices from
  !!        any subdivided neighbour (collectEdgeVerts) so the fan threads through
  !!        every vertex the neighbour introduced along the shared edge.
  !!   5. Ray cast against the resulting watertight triangle mesh (Woop et al. 2013).
  !!
  !! Patch adjacency is precomputed at init by matching shared boundary control points.
  !!
  !! Edges of each patch (indices for adj array):
  !!   1 = u=0 boundary: row j=1, varying in v (C00 to C01)
  !!   2 = u=1 boundary: row j=4, varying in v (C10 to C11)
  !!   3 = v=0 boundary: col k=1, varying in u (C00 to C10)
  !!   4 = v=1 boundary: col k=4, varying in u (C01 to C11)
  !!
  !! Sample input:
  !!   vol { type bezierVolume; id 1; numPatches 6; ctrlPts (x y z ...); }
  !!
  !! See misclassClerk_class (Tallies/TallyClerks) for a halfspace
  !! misclassification diagnostic against a reference region.
  !!
  !! Private Members:
  !!   ctrlPts    -> Control points of every patch, (numPatches, 4, 4, 3)
  !!   numPatches -> Number of patches
  !!   allPtsFlat -> All control points flattened to one list, for the global GJK
  !!     rejection test (Step 2 of the halfspace algorithm above)
  !!   nAllPts    -> Number of points in allPtsFlat
  !!   aabb       -> Axis-aligned bounding box of the whole volume, (xMin,yMin,zMin,xMax,yMax,zMax)
  !!   adj        -> adj(i,e) = index of the patch sharing edge e of patch i; -1 if none
  !!   adjEdge    -> adjEdge(i,e) = edge index on that adjacent patch; -1 if none
  !!   weights    -> Rational (NURBS) weight of every control point, (numPatches, 4, 4)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: bezierVolume
    private
    real(defReal), dimension(:,:,:,:), allocatable :: ctrlPts      ! (numPatches, 4, 4, 3)
    integer(shortInt)                              :: numPatches = 0
    real(defReal), dimension(:,:), allocatable     :: allPtsFlat   ! (nAllPts, 3) for global GJK
    integer(shortInt)                              :: nAllPts    = 0
    real(defReal), dimension(6)                    :: aabb       = ZERO
    ! adj(i,e)     = index of patch sharing edge e of patch i; -1 if none
    ! adjEdge(i,e) = edge index on that adjacent patch; -1 if none
    integer(shortInt), dimension(:,:), allocatable :: adj
    integer(shortInt), dimension(:,:), allocatable :: adjEdge
    real(defReal), dimension(:,:,:), allocatable   :: weights      ! (numPatches, 4, 4)
  contains
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
    procedure :: buildAdjacency
    procedure :: halfspace
    procedure :: inPatchAABB
    procedure :: rayCast
    procedure :: rayTriangle
  end type bezierVolume


contains

  ! ---------------------------------------------------------------------------
  ! Surface interface routines
  ! ---------------------------------------------------------------------------

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
  !!   fatalError if the number of ctrlPts entries does not match numPatches * 16 * 3
  !!   fatalError if the number of weights entries (when present) does not match
  !!     numPatches * 16
  !!   fatalError (via buildAdjacency) if the input patches do not form a
  !!     watertight closed volume
  !!
  subroutine init(self, dict)
    class(bezierVolume), intent(inout)       :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, n, m, i, j, k, l
    real(defReal), dimension(:), allocatable :: ctrlPtsList, weightsList
    character(100), parameter :: Here = 'init (bezierVolume_class.f90)'

    call dict % get(id,              'id')
    call dict % get(ctrlPtsList,     'ctrlPts')
    call dict % get(self % numPatches, 'numPatches')

    if (id < 1) call fatalError(Here, 'Invalid surface id. Must be >= 1')

    n = size(ctrlPtsList)
    if (n /= self % numPatches * 16 * 3) then
      call fatalError(Here, 'Control points inconsistent with numPatches. '// &
                            'Expected: '//numToChar(self % numPatches * 16 * 3)// &
                            ' Got: '//numToChar(n))
    end if

    call self % setID(id)

    ! Load control points into (numPatches, 4, 4, 3) array
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

    ! Load weights (optional; defaults to 1.0 for polynomial Bezier)
    allocate(self % weights(self % numPatches, 4, 4))
    if (dict % isPresent('weights')) then
      call dict % get(weightsList, 'weights')
      n = size(weightsList)
      if (n /= self % numPatches * 16) then
        call fatalError(Here, 'Weights inconsistent with numPatches. '// &
                              'Expected: '//numToChar(self % numPatches * 16)// &
                              ' Got: '//numToChar(n))
      end if
      m = 1
      do i = 1, self % numPatches
        do j = 1, 4
          do k = 1, 4
            self % weights(i, j, k) = weightsList(m)
            m = m + 1
          end do
        end do
      end do
    else
      self % weights = ONE
    end if

    self % aabb = self % boundingBox()

    ! Flat array of all control points for the global GJK test
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

    call self % buildAdjacency()

  end subroutine init

  !!
  !! Return axis-aligned bounding box of the volume
  !!
  !! See surface_inter for more details
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
  !! Permanent no-op stub: this class overrides halfspace directly instead of
  !! using the base surface class's default evaluate/going-driven halfspace,
  !! so evaluate is never actually called and has no meaningful value to return.
  !!
  !! See surface_inter for more details
  !!
  pure function evaluate(self, r) result(c)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    c = ZERO
  end function evaluate

  !!
  !! Permanent no-op stub: this class overrides halfspace directly (see
  !! evaluate above), so distance-to-surface is never queried through here.
  !!
  !! See surface_inter for more details
  !!
  pure function distance(self, r, u) result(d)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: d
    d = INF
  end function distance

  !!
  !! Permanent no-op stub: this class overrides halfspace directly (see
  !! evaluate above), so this generic evaluate-sign convention is never used.
  !!
  !! See surface_inter for more details
  !!
  pure function going(self, r, u) result(halfspace)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: halfspace
    halfspace = .false.
  end function going

  !!
  !! Return to uninitialised state
  !!
  !! See surface_inter for more details
  !!
  elemental subroutine kill(self)
    class(bezierVolume), intent(inout) :: self

    call kill_super(self)

    if (allocated(self % ctrlPts))    deallocate(self % ctrlPts)
    if (allocated(self % weights))    deallocate(self % weights)
    if (allocated(self % allPtsFlat)) deallocate(self % allPtsFlat)
    if (allocated(self % adj))        deallocate(self % adj)
    if (allocated(self % adjEdge))    deallocate(self % adjEdge)

    self % numPatches = 0
    self % nAllPts    = 0
    self % aabb       = ZERO

  end subroutine kill

  ! ---------------------------------------------------------------------------
  ! Adjacency construction
  ! ---------------------------------------------------------------------------

  !!
  !! Checks whether all 4 corners of a patch edge have collapsed to one point
  !! (e.g. a polar/fan apex), to within SURF_TOL.
  !!
  pure function isDegenerateEdge(edge) result(deg)
    real(defReal), dimension(4,3), intent(in) :: edge
    logical(defBool) :: deg
    integer(shortInt) :: k

    deg = .true.
    do k = 2, 4
      if (any(abs(edge(k,:) - edge(1,:)) > SURF_TOL)) deg = .false.
    end do

  end function isDegenerateEdge

  !!
  !! Build the patch adjacency map self%adj(numPatches, 4).
  !! Two patches share an edge if all 4 boundary control points match
  !! (forward or reverse order) to within SURF_TOL.
  !!
  !! Degenerate edges (all 4 corners collapsed to one point, e.g. a polar/fan
  !! apex) are exempted from matching: such a point can be shared by more
  !! than 2 patches at once, but this loop only ever records one neighbour
  !! per edge, so matching a degenerate edge would leave an inconsistent,
  !! asymmetric adjacency graph. A degenerate edge has no interior length and
  !! therefore no T-junction vertices to reconcile, so leaving it unmatched
  !! (adj/adjEdge = -1) is safe -- buildWatertightMesh's own AABB-fallback
  !! logic already treats that as "no adjacent patch".
  !!
  !! Errors:
  !!   fatalError (via checkWatertight) if the input patches do not form a
  !!     watertight closed volume
  !!
  subroutine buildAdjacency(self)
    class(bezierVolume), intent(inout) :: self
    integer(shortInt) :: i, j, ei, ej, k
    real(defReal), dimension(4,3) :: edgeI, edgeJ
    logical(defBool) :: fwd, rev

    allocate(self % adj(self % numPatches, 4))
    allocate(self % adjEdge(self % numPatches, 4))
    self % adj     = -1
    self % adjEdge = -1

    do i = 1, self % numPatches
      do ei = 1, 4
        edgeI = patchEdge(self % ctrlPts(i,:,:,:), ei)
        if (isDegenerateEdge(edgeI)) cycle

        do j = i + 1, self % numPatches
          do ej = 1, 4
            edgeJ = patchEdge(self % ctrlPts(j,:,:,:), ej)
            fwd = .true.
            rev = .true.
            do k = 1, 4
              if (any(abs(edgeI(k,:) - edgeJ(k,:))   > SURF_TOL)) fwd = .false.
              if (any(abs(edgeI(k,:) - edgeJ(5-k,:)) > SURF_TOL)) rev = .false.
            end do
            if (fwd .or. rev) then
              self % adj(i, ei)     = j;  self % adjEdge(i, ei) = ej
              self % adj(j, ej)     = i;  self % adjEdge(j, ej) = ei
            end if
          end do
        end do
      end do
    end do

    call checkWatertight(self)

  end subroutine buildAdjacency

  !!
  !! Check that every non-degenerate patch edge found a matching neighbour in
  !! buildAdjacency's matching loop above. A dangling edge means the input
  !! patches do not form a closed volume, and would otherwise surface only
  !! much later as an unexplained halfspace misclassification -- so it is
  !! reported here, at construction time.
  !!
  !! Errors:
  !!   fatalError if any non-degenerate edge has no matching neighbour
  !!
  subroutine checkWatertight(self)
    class(bezierVolume), intent(in) :: self
    integer(shortInt) :: i, ei
    real(defReal), dimension(4,3) :: edgeI
    character(100), parameter :: Here = 'checkWatertight (bezierVolume_class.f90)'

    do i = 1, self % numPatches
      do ei = 1, 4
        if (self % adj(i, ei) /= -1) cycle
        edgeI = patchEdge(self % ctrlPts(i,:,:,:), ei)
        if (isDegenerateEdge(edgeI)) cycle
        call fatalError(Here, 'patch '//numToChar(i)//' edge '//numToChar(ei)// &
                               ' has no matching neighbour: geometry is not watertight')
      end do
    end do

  end subroutine checkWatertight

  !!
  !! Extract the 4 control points of edge e from a patch:
  !!   e=1: u=0 boundary, row j=1, varying in v (C00 -> C01)
  !!   e=2: u=1 boundary, row j=4, varying in v (C10 -> C11)
  !!   e=3: v=0 boundary, col k=1, varying in u (C00 -> C10)
  !!   e=4: v=1 boundary, col k=4, varying in u (C01 -> C11)
  !!
  pure function patchEdge(patch, e) result(edge)
    real(defReal), dimension(4,4,3), intent(in) :: patch
    integer(shortInt), intent(in)               :: e
    real(defReal), dimension(4,3)               :: edge
    integer(shortInt) :: k

    select case(e)
      case(1);  do k = 1, 4;  edge(k,:) = patch(1,k,:);  end do
      case(2);  do k = 1, 4;  edge(k,:) = patch(4,k,:);  end do
      case(3);  do k = 1, 4;  edge(k,:) = patch(k,1,:);  end do
      case(4);  do k = 1, 4;  edge(k,:) = patch(k,4,:);  end do
      case default;  edge = ZERO
    end select

  end function patchEdge

  ! ---------------------------------------------------------------------------
  ! Halfspace pipeline
  ! ---------------------------------------------------------------------------

  !!
  !! Determine halfspace for a particle position.
  !!
  !! Result: .true. = outside (positive halfspace)
  !!
  function halfspace(self, r, u) result(hs)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: hs

    type(subPatches) :: subPatch
    real(defReal), dimension(:,:,:),   allocatable :: curWeights
    ! Per ORIGINAL patch: whether it was split at all this call
    logical(defBool), dimension(self % numPatches) :: wasSubdivided
    real(defReal), dimension(:,:,:), allocatable :: tris
    integer(shortInt) :: nTri

    ! --- Step 1: AABB rejection ---
    if (r(1) < self % aabb(1) - SURF_TOL .or. r(1) > self % aabb(4) + SURF_TOL .or. &
        r(2) < self % aabb(2) - SURF_TOL .or. r(2) > self % aabb(5) + SURF_TOL .or. &
        r(3) < self % aabb(3) - SURF_TOL .or. r(3) > self % aabb(6) + SURF_TOL) then
      hs = .true.
      return
    end if

    ! --- Step 2: Global convex hull rejection ---
    if (.not. pointInConvexHull(self % allPtsFlat, self % nAllPts, r)) then
      hs = .true.
      return
    end if

    ! --- Step 3: Selective subdivision ---
    call subdivide(self, r, subPatch, curWeights, wasSubdivided)

    ! --- Step 4: Build watertight triangle mesh with T-junction patching ---
    call buildWatertightMesh(self, subPatch, wasSubdivided, tris, nTri)

    ! --- Step 5: Ray cast against the watertight triangle mesh ---
    hs = .not. self % rayCast(tris, nTri, r)

    deallocate(curWeights, tris)

  end function halfspace

  !!
  !! Adaptively subdivide the patches that may contain query point r, via AABB
  !! + GJK convex-hull containment on each candidate sub-patch, up to
  !! MAX_SUBDIVISIONS levels deep. Patches whose bounding volume cannot
  !! contain r are left unsplit.
  !!
  !! Args:
  !!   self [in] -> the bezierVolume
  !!   r [in] -> query point
  !!   subPatch [out] -> the resulting sub-patch list
  !!   curWeights [out] -> sub-patch weights, curWeights(1:subPatch%n,4,4)
  !!   wasSubdivided [out] -> per original patch, whether it was split at all
  !!
  subroutine subdivide(self, r, subPatch, curWeights, wasSubdivided)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    type(subPatches), intent(out) :: subPatch
    real(defReal), dimension(:,:,:),   allocatable, intent(out) :: curWeights
    logical(defBool), dimension(self % numPatches), intent(out) :: wasSubdivided

    real(defReal), dimension(:,:,:,:), allocatable :: newPts
    real(defReal), dimension(:,:,:),   allocatable :: newWeights
    integer(shortInt), dimension(:), allocatable   :: newOrig
    logical(defBool), dimension(:), allocatable    :: needsSubdiv
    real(defReal), dimension(:,:), allocatable     :: newUV
    ! Original patches that have passed the STRICT hull test at least once this
    ! call. Used to scope the AABB-only fallback below to patches actually
    ! adjacent to the region of interest, instead of scanning every patch in the
    ! model by loose bounding-box coincidence.
    logical(defBool), dimension(self % numPatches) :: everHullSelected
    ! Candidate original-patch index for the AABB-fallback eligibility check
    integer(shortInt) :: candidatePatch
    real(defReal), dimension(16, 3) :: patchPts
    real(defReal), dimension(4,4,3) :: Q00, Q01, Q10, Q11
    real(defReal), dimension(4,4)   :: Q00w, Q01w, Q10w, Q11w
    integer(shortInt) :: subdivision, i, j, k, m, newN, nNeed
    real(defReal) :: u0, u1, v0, v1, um, vm

    subPatch % n = self % numPatches
    allocate(subPatch % ctrlPts(subPatch % n, 4, 4, 3))
    allocate(curWeights(subPatch % n, 4, 4))
    allocate(subPatch % orig(subPatch % n))
    allocate(subPatch % uvRange(subPatch % n, 4))
    subPatch % ctrlPts = self % ctrlPts
    curWeights   = self % weights
    do i = 1, subPatch % n
      subPatch % orig(i) = i
      subPatch % uvRange(i, :) = (/ ZERO, ONE, ZERO, ONE /)
    end do
    wasSubdivided = .false.
    everHullSelected = .false.

    do subdivision = 1, MAX_SUBDIVISIONS

      ! Identify which current patches need subdivision (AABB + GJK on 16 control points)
      allocate(needsSubdiv(subPatch % n))
      needsSubdiv = .false.
      nNeed = 0

      do i = 1, subPatch % n
        if (.not. self % inPatchAABB(subPatch % ctrlPts(i,:,:,:), r)) cycle
        m = 0
        do j = 1, 4
          do k = 1, 4
            m = m + 1
            patchPts(m,:) = subPatch % ctrlPts(i, j, k, :)
          end do
        end do
        if (pointInConvexHull(patchPts, 16, r)) then
          needsSubdiv(i) = .true.
          nNeed = nNeed + 1
          everHullSelected(subPatch % orig(i)) = .true.
        end if
      end do

      if (nNeed == 0) then
        ! GJK found no patch (typically floating-point precision on an
        ! already-tight, deeply-subdivided hull): fall back to AABB-only
        ! selection, scoped via isFallbackEligible (see its doc comment).
        do i = 1, subPatch % n
          if (.not. self % inPatchAABB(subPatch % ctrlPts(i,:,:,:), r)) cycle
          candidatePatch = subPatch % orig(i)
          if (isFallbackEligible(candidatePatch, everHullSelected, self % adj)) then
            needsSubdiv(i) = .true.
            nNeed = nNeed + 1
          end if
        end do
        if (nNeed == 0) then
          deallocate(needsSubdiv)
          exit
        end if
      end if

      ! Build replacement list: keep unsplit patches, expand split ones into 4 each
      newN = subPatch % n + 3 * nNeed
      allocate(newPts(newN, 4, 4, 3))
      allocate(newWeights(newN, 4, 4))
      allocate(newOrig(newN))
      allocate(newUV(newN, 4))
      j = 0
      do i = 1, subPatch % n
        if (needsSubdiv(i)) then
          wasSubdivided(subPatch % orig(i)) = .true.
          call subdividePatch(subPatch % ctrlPts(i,:,:,:), curWeights(i,:,:), HALF, HALF, &
                              Q00, Q01, Q10, Q11, Q00w, Q01w, Q10w, Q11w)
          u0 = subPatch % uvRange(i,1);  u1 = subPatch % uvRange(i,2)
          v0 = subPatch % uvRange(i,3);  v1 = subPatch % uvRange(i,4)
          um = HALF*(u0+u1);  vm = HALF*(v0+v1)
          ! Q00: u∈[u0,um], v∈[v0,vm]
          j = j + 1;  newPts(j,:,:,:) = Q00;  newWeights(j,:,:) = Q00w
          newOrig(j) = subPatch % orig(i);  newUV(j,:) = (/ u0, um, v0, vm /)
          ! Q01: u∈[u0,um], v∈[vm,v1]
          j = j + 1;  newPts(j,:,:,:) = Q01;  newWeights(j,:,:) = Q01w
          newOrig(j) = subPatch % orig(i);  newUV(j,:) = (/ u0, um, vm, v1 /)
          ! Q10: u∈[um,u1], v∈[v0,vm]
          j = j + 1;  newPts(j,:,:,:) = Q10;  newWeights(j,:,:) = Q10w
          newOrig(j) = subPatch % orig(i);  newUV(j,:) = (/ um, u1, v0, vm /)
          ! Q11: u∈[um,u1], v∈[vm,v1]
          j = j + 1;  newPts(j,:,:,:) = Q11;  newWeights(j,:,:) = Q11w
          newOrig(j) = subPatch % orig(i);  newUV(j,:) = (/ um, u1, vm, v1 /)
        else
          j = j + 1
          newPts(j,:,:,:)  = subPatch % ctrlPts(i,:,:,:)
          newWeights(j,:,:) = curWeights(i,:,:)
          newOrig(j)        = subPatch % orig(i)
          newUV(j,:)        = subPatch % uvRange(i,:)
        end if
      end do

      deallocate(needsSubdiv)
      call move_alloc(newPts,     subPatch % ctrlPts)
      call move_alloc(newWeights, curWeights)
      call move_alloc(newOrig,    subPatch % orig)
      call move_alloc(newUV,      subPatch % uvRange)
      subPatch % n = newN

    end do

  end subroutine subdivide

  !!
  !! Build a watertight triangle mesh over the whole bezierVolume from the
  !! (possibly non-uniformly subdivided) sub-patch list produced by
  !! subdivide, patching T-junctions between sub-patches at
  !! different subdivision depths -- both within one original patch, and
  !! across original patches where one side was subdivided and its
  !! neighbour was not.
  !!
  !! Dispatches each original patch to triSubPatch or
  !! triUnsubPatch, which document the two triangulation
  !! strategies in full.
  !!
  !! Triangle orientation (winding unused — we just count ray-cast hits):
  !!   T1 = (C00, C01, C11)   covers the u=0 and v=1 boundaries
  !!   T2 = (C00, C11, C10)   covers the u=1 and v=0 boundaries
  !!
  !! Args:
  !!   self [in] -> the bezierVolume
  !!   subPatch [in] -> the current sub-patch list, from subdivide
  !!   wasSubdivided [in] -> per original patch, whether it was split at all
  !!   tris [out] -> triangle mesh, tris(1:nTri,3,3)
  !!   nTri [out] -> number of valid triangles in tris
  !!
  subroutine buildWatertightMesh(self, subPatch, wasSubdivided, tris, nTri)
    class(bezierVolume), intent(in) :: self
    type(subPatches), intent(in)    :: subPatch
    logical(defBool), dimension(self % numPatches), intent(in) :: wasSubdivided
    real(defReal), dimension(:,:,:), allocatable, intent(out)  :: tris
    integer(shortInt), intent(out) :: nTri

    integer(shortInt) :: P, maxTri

    maxTri = MAX_TJ_BUF * subPatch % n + MAX_TJ_BUF * self % numPatches
    allocate(tris(maxTri, 3, 3))
    nTri = 0

    do P = 1, self % numPatches
      if (wasSubdivided(P)) then
        call triSubPatch(self, P, subPatch, wasSubdivided, tris, nTri)
      else
        call triUnsubPatch(self, P, subPatch, wasSubdivided, tris, nTri)
      end if
    end do

  end subroutine buildWatertightMesh

  !!
  !! Polygon-fan triangulation for a subdivided original patch P: fan each of
  !! its sub-patches, closing two kinds of T-junction gap.
  !!
  !! Within-patch: selective subdivision may leave adjacent sub-patches of the
  !! same original patch at different depths (collectIntraPatchEdgeVerts).
  !!
  !! Cross-patch, along P's true outer boundary: when P's edge there is
  !! shared with an ALSO-subdivided original patch, neither side's
  !! independent quadtree refinement is otherwise synchronised with the
  !! other, so the neighbour's boundary vertices are pulled in
  !! (collectEdgeVerts) and merged into the same buffer used for the fan
  !! (mergeEdgeVerts), giving both sides' triangulations an identical vertex
  !! set. Skipped for edges internal to P (shared with another sub-patch of
  !! the same original patch), since those are already handled by the
  !! within-patch case above.
  !!
  !! Args:
  !!   self [in] -> the bezierVolume
  !!   P [in] -> the original patch index being triangulated
  !!   subPatch [in] -> the current sub-patch list, from subdivide
  !!   wasSubdivided [in] -> per original patch, whether it was split at all
  !!   tris [inout], nTri [inout] -> triangle mesh being built
  !!
  subroutine triSubPatch(self, P, subPatch, wasSubdivided, tris, nTri)
    class(bezierVolume), intent(in) :: self
    integer(shortInt), intent(in)   :: P
    type(subPatches), intent(in)    :: subPatch
    logical(defBool), dimension(self % numPatches), intent(in) :: wasSubdivided
    real(defReal), dimension(:,:,:), intent(inout) :: tris
    integer(shortInt), intent(inout) :: nTri

    integer(shortInt) :: i
    real(defReal), dimension(3) :: C00, C01, C10, C11
    real(defReal), dimension(MAX_TJ_BUF, 3) :: e1buf, e2buf, e3buf, e4buf
    integer(shortInt) :: n1, n2, n3, n4

    do i = 1, subPatch % n
      if (subPatch % orig(i) /= P) cycle
      C00 = subPatch % ctrlPts(i, 1, 1, :)
      C01 = subPatch % ctrlPts(i, 1, 4, :)
      C10 = subPatch % ctrlPts(i, 4, 1, :)
      C11 = subPatch % ctrlPts(i, 4, 4, :)

      call collectIntraPatchEdgeVerts(subPatch, i, 1, C00, C01, e1buf, n1)
      call collectIntraPatchEdgeVerts(subPatch, i, 2, C10, C11, e2buf, n2)
      call collectIntraPatchEdgeVerts(subPatch, i, 3, C00, C10, e3buf, n3)
      call collectIntraPatchEdgeVerts(subPatch, i, 4, C01, C11, e4buf, n4)

      if (edgeOnOuterBoundary(subPatch % uvRange, i, 1)) then
        call reconcileEdge(self, P, 1, C00, C01, subPatch, wasSubdivided, e1buf, n1)
      end if
      if (edgeOnOuterBoundary(subPatch % uvRange, i, 2)) then
        call reconcileEdge(self, P, 2, C10, C11, subPatch, wasSubdivided, e2buf, n2)
      end if
      if (edgeOnOuterBoundary(subPatch % uvRange, i, 3)) then
        call reconcileEdge(self, P, 3, C00, C10, subPatch, wasSubdivided, e3buf, n3)
      end if
      if (edgeOnOuterBoundary(subPatch % uvRange, i, 4)) then
        call reconcileEdge(self, P, 4, C01, C11, subPatch, wasSubdivided, e4buf, n4)
      end if

      ! T1: fan over polygon C00 → [e1 interior] → C01 → [e4 interior] → C11
      call fanTriangulate(e1buf, n1, e4buf, n4, tris, nTri)

      ! T2: fan over polygon C00 → [e3 interior] → C10 → [e2 interior] → C11
      call fanTriangulate(e3buf, n3, e2buf, n2, tris, nTri)
    end do

  end subroutine triSubPatch

  !!
  !! Reconcile sub-patch i's outer-boundary edge edgeIdx (spanning
  !! cornerA-cornerB) against an ALSO-subdivided original patch adjacent to P
  !! along that edge, merging the neighbour's boundary vertices into
  !! buf/nBuf. No-op if there is no adjacent patch, or if it was not itself
  !! subdivided. See triSubPatch's own doc comment for why
  !! this reconciliation is needed.
  !!
  !! Args:
  !!   self [in] -> the bezierVolume
  !!   P [in] -> the original patch index
  !!   edgeIdx [in] -> which of P's 4 edges (1-4) this is
  !!   cornerA [in], cornerB [in] -> the two corner points of edge edgeIdx
  !!   subPatch [in] -> the current sub-patch list, from subdivide
  !!   wasSubdivided [in] -> per original patch, whether it was split at all
  !!   buf [inout], nBuf [inout] -> this sub-patch's edge-vertex buffer to
  !!     merge the neighbour's vertices into
  !!
  subroutine reconcileEdge(self, P, edgeIdx, cornerA, cornerB, subPatch, wasSubdivided, buf, nBuf)
    class(bezierVolume), intent(in)         :: self
    integer(shortInt), intent(in)           :: P, edgeIdx
    real(defReal), dimension(3), intent(in) :: cornerA, cornerB
    type(subPatches), intent(in)            :: subPatch
    logical(defBool), dimension(self % numPatches), intent(in) :: wasSubdivided
    real(defReal), dimension(:,:), intent(inout)  :: buf
    integer(shortInt), intent(inout)              :: nBuf

    real(defReal), dimension(MAX_TJ_BUF, 3) :: xBuf
    integer(shortInt) :: nX, adjPatch, adjEdgeIdx

    adjPatch = self % adj(P, edgeIdx)
    if (adjPatch <= 0) return
    if (.not. wasSubdivided(adjPatch)) return

    adjEdgeIdx = self % adjEdge(P, edgeIdx)
    call collectEdgeVerts(subPatch, adjPatch, adjEdgeIdx, cornerA, cornerB, xBuf, nX)
    call mergeEdgeVerts(buf, nBuf, xBuf, nX, cornerA, cornerB)

  end subroutine reconcileEdge

  !!
  !! Polygon-fan triangulation for an unsubdivided original patch Q, with
  !! full T-junction matching: for each boundary edge of Q adjacent to a
  !! subdivided patch, collect ALL sub-patch corners lying on that shared
  !! edge (collectEdgeVerts), sorted along the edge, and fan-triangulate from
  !! C00 through them -- closing multi-level T-junction gaps rather than
  !! just a single midpoint.
  !!
  !! Args: as triSubPatch.
  !!
  subroutine triUnsubPatch(self, P, subPatch, wasSubdivided, tris, nTri)
    class(bezierVolume), intent(in) :: self
    integer(shortInt), intent(in)   :: P
    type(subPatches), intent(in)    :: subPatch
    logical(defBool), dimension(self % numPatches), intent(in) :: wasSubdivided
    real(defReal), dimension(:,:,:), intent(inout) :: tris
    integer(shortInt), intent(inout) :: nTri

    real(defReal), dimension(3) :: C00, C01, C10, C11
    real(defReal), dimension(MAX_TJ_BUF, 3) :: e1buf, e2buf, e3buf, e4buf
    integer(shortInt) :: n1, n2, n3, n4

    C00 = self % ctrlPts(P, 1, 1, :)
    C01 = self % ctrlPts(P, 1, 4, :)
    C10 = self % ctrlPts(P, 4, 1, :)
    C11 = self % ctrlPts(P, 4, 4, :)

    call collectOuters(self, P, 1, C00, C01, subPatch, wasSubdivided, e1buf, n1)
    call collectOuters(self, P, 2, C10, C11, subPatch, wasSubdivided, e2buf, n2)
    call collectOuters(self, P, 3, C00, C10, subPatch, wasSubdivided, e3buf, n3)
    call collectOuters(self, P, 4, C01, C11, subPatch, wasSubdivided, e4buf, n4)

    ! T1: fan over polygon C00 → [e1 interior] → C01 → [e4 interior] → C11
    call fanTriangulate(e1buf, n1, e4buf, n4, tris, nTri)

    ! T2: fan over polygon C00 → [e3 interior] → C10 → [e2 interior] → C11
    call fanTriangulate(e3buf, n3, e2buf, n2, tris, nTri)

  end subroutine triUnsubPatch

  !!
  !! Collect the T-junction vertices along P's edge edgeIdx (spanning
  !! cornerA-cornerB), when P itself was NOT subdivided: if the adjacent
  !! original patch along that edge WAS subdivided, collect all its boundary
  !! vertices via collectEdgeVerts (sorted along the edge); otherwise the
  !! edge has no T-junctions and buf is just its two corners.
  !!
  !! Args:
  !!   self [in] -> the bezierVolume
  !!   P [in] -> the original patch index
  !!   edgeIdx [in] -> which of P's 4 edges (1-4) this is
  !!   cornerA [in], cornerB [in] -> the two corner points of edge edgeIdx
  !!   subPatch [in] -> the current sub-patch list, from subdivide
  !!   wasSubdivided [in] -> per original patch, whether it was split at all
  !!   buf [out] -> collected edge-vertex buffer
  !!   nBuf [out] -> number of valid entries in buf
  !!
  subroutine collectOuters(self, P, edgeIdx, cornerA, cornerB, subPatch, wasSubdivided, buf, nBuf)
    class(bezierVolume), intent(in)         :: self
    integer(shortInt), intent(in)           :: P, edgeIdx
    real(defReal), dimension(3), intent(in) :: cornerA, cornerB
    type(subPatches), intent(in)            :: subPatch
    logical(defBool), dimension(self % numPatches), intent(in) :: wasSubdivided
    real(defReal), dimension(:,:), intent(out) :: buf
    integer(shortInt), intent(out)             :: nBuf

    integer(shortInt) :: adjPatch, adjEdgeIdx

    adjPatch = self % adj(P, edgeIdx)
    if (adjPatch > 0 .and. wasSubdivided(adjPatch)) then
      adjEdgeIdx = self % adjEdge(P, edgeIdx)
      call collectEdgeVerts(subPatch, adjPatch, adjEdgeIdx, cornerA, cornerB, buf, nBuf)
    else
      nBuf = 2;  buf(1,:) = cornerA;  buf(2,:) = cornerB
    end if

  end subroutine collectOuters

  !!
  !! AABB pre-filter for a single patch's 16 control points.
  !!
  pure function inPatchAABB(self, patch, r) result(inside)
    class(bezierVolume), intent(in)             :: self
    real(defReal), dimension(4,4,3), intent(in) :: patch
    real(defReal), dimension(3), intent(in)     :: r
    logical(defBool)                            :: inside
    real(defReal), dimension(3)                 :: minPt, maxPt
    integer(shortInt)                           :: j, k

    minPt =  INF;  maxPt = -INF
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
  !! Checks whether an original patch is eligible for the AABB-only subdivision fallback
  !!
  !! Used by subdivide's AABB-only fallback, which triggers when
  !! the strict GJK hull test selects no patch at all (typically
  !! floating-point precision on an already-tight, deeply-subdivided hull).
  !! Scoping the fallback matters because an unscoped AABB scan pulls in
  !! totally unrelated patches purely by loose bounding-box coincidence: a
  !! patch can be selected despite its own hull never once containing the
  !! query point, then only ever refined one level deep while its genuine
  !! neighbour keeps refining much further -- a mismatch that opens a real
  !! gap between them. Eligible patches are restricted to ones that have
  !! themselves passed the strict hull test at least once this call, or are
  !! adjacent to one that has. If no patch has passed the strict test at all
  !! yet, there is no adjacency information to scope by, so every patch is
  !! eligible (matches the original unrestricted-scan behaviour for that
  !! narrow case).
  !!
  !! Args:
  !!   candidatePatch [in] -> index of the original patch being considered for
  !!     the fallback
  !!   everHullSelected [in] -> per-original-patch flag, true if that patch has
  !!     passed the strict hull test at least once this halfspace() call
  !!   adj [in] -> patch adjacency map, adj(i,e) = index of the patch sharing
  !!     edge e of patch i, or -1 if none
  !!
  !! Result:
  !!   .true. if no original patch has EVER passed the strict hull test this
  !!   call (no scoping information yet), if candidatePatch itself has, or if
  !!   one of candidatePatch's 4 neighbours has; .false. otherwise.
  !!
  pure function isFallbackEligible(candidatePatch, everHullSelected, adj) result(eligible)
    integer(shortInt), intent(in)                 :: candidatePatch
    logical(defBool), dimension(:), intent(in)    :: everHullSelected
    integer(shortInt), dimension(:,:), intent(in) :: adj
    logical(defBool)                              :: eligible
    integer(shortInt) :: e, adjPatch

    if (.not. any(everHullSelected)) then
      eligible = .true.
    else if (everHullSelected(candidatePatch)) then
      eligible = .true.
    else
      eligible = .false.
      do e = 1, 4
        adjPatch = adj(candidatePatch, e)
        if (adjPatch > 0) then
          if (everHullSelected(adjPatch)) eligible = .true.
        end if
      end do
    end if

  end function isFallbackEligible

  !!
  !! Append pt to buf(1:nBuf) unless a point within tol is already present.
  !!
  !! Args:
  !!   buf [inout] -> point buffer
  !!   nBuf [inout] -> number of valid entries in buf; incremented on insert
  !!   pt [in] -> candidate point
  !!   tol [in] -> coincidence tolerance
  !!
  !! Errors:
  !!   fatalError if buf is already at capacity when a new point would be inserted
  !!
  subroutine insertUnique(buf, nBuf, pt, tol)
    real(defReal), dimension(:,:), intent(inout) :: buf
    integer(shortInt), intent(inout)             :: nBuf
    real(defReal), dimension(3), intent(in)      :: pt
    real(defReal), intent(in)                    :: tol
    integer(shortInt) :: si
    logical(defBool)  :: isDup
    character(100), parameter :: Here = 'insertUnique (bezierVolume_class.f90)'

    isDup = .false.
    do si = 1, nBuf
      if (norm2(buf(si,:) - pt) < tol) then
        isDup = .true.
        exit
      end if
    end do
    if (.not. isDup) then
      if (nBuf >= size(buf,1)) call fatalError(Here, 'buffer full ('//numToChar(size(buf,1))//')')
      nBuf = nBuf + 1
      buf(nBuf,:) = pt
    end if

  end subroutine insertUnique

  !!
  !! Sort buf(1:nBuf) in place by ascending key(1:nBuf) (insertion sort;
  !! buf holds at most a few tens of entries, so O(n^2) is not a concern).
  !!
  !! Args:
  !!   buf [inout] -> point buffer, reordered in place
  !!   nBuf [in] -> number of valid entries in buf
  !!   key [inout] -> sort key, reordered alongside buf
  !!
  subroutine sortByKey(buf, nBuf, key)
    real(defReal), dimension(:,:), intent(inout) :: buf
    integer(shortInt), intent(in)                :: nBuf
    real(defReal), dimension(:), intent(inout)   :: key
    integer(shortInt) :: i, si
    real(defReal), dimension(3) :: tmpV
    real(defReal) :: tmpD

    do si = 2, nBuf
      tmpV = buf(si,:);  tmpD = key(si)
      i = si - 1
      do while (i >= 1 .and. key(i) > tmpD)
        buf(i+1,:) = buf(i,:);  key(i+1) = key(i)
        i = i - 1
      end do
      buf(i+1,:) = tmpV;  key(i+1) = tmpD
    end do

  end subroutine sortByKey

  !!
  !! Collect T-junction vertices on one edge of sub-patch patchIdx from finer
  !! sub-patches of the same original patch (within-patch T-junction detection).
  !!
  !! edgeIdx follows the patch edge numbering:
  !!   1 = u=uMin (left),  2 = u=uMax (right)
  !!   3 = v=vMin (bottom), 4 = v=vMax (top)
  !!
  !! ptA and ptB are the 3D end-points of patchIdx's edge in the direction the
  !! polygon fan expects them (matches the T1/T2 polygon winding in halfspace).
  !!
  !! Returns buf sorted ptA→ptB with all T-junction vertices inserted.
  !! Always contains at least ptA and ptB.
  !!
  subroutine collectIntraPatchEdgeVerts(subPatch, patchIdx, edgeIdx, ptA, ptB, buf, nBuf)
    type(subPatches), intent(in)            :: subPatch
    integer(shortInt), intent(in)           :: patchIdx, edgeIdx
    real(defReal), dimension(3), intent(in) :: ptA, ptB
    real(defReal), dimension(:,:), intent(out) :: buf
    integer(shortInt), intent(out)             :: nBuf

    integer(shortInt) :: i, P, nBufBefore
    real(defReal), dimension(3)       :: pt1, pt2
    real(defReal), dimension(size(buf,1)) :: paramArr
    real(defReal) :: boundaryVal, pMin, pMax, p1, p2
    logical(defBool) :: onEdge

    P    = subPatch % orig(patchIdx)
    nBuf = 0
    ! 2026-08-22 fix: sort by true Bezier parameter, not chord projection --
    ! see the matching note in collectEdgeVerts. No direction ambiguity here:
    ! Si and Sj are sub-patches of the SAME original patch P, so increasing u/v
    ! always runs ptA->ptB (subdivision never reverses an axis' orientation).

    ! UV boundary value and range for this edge of Si
    select case (edgeIdx)
      case(1);  boundaryVal = subPatch % uvRange(patchIdx,1)
                pMin = subPatch % uvRange(patchIdx,3);  pMax = subPatch % uvRange(patchIdx,4)
      case(2);  boundaryVal = subPatch % uvRange(patchIdx,2)
                pMin = subPatch % uvRange(patchIdx,3);  pMax = subPatch % uvRange(patchIdx,4)
      case(3);  boundaryVal = subPatch % uvRange(patchIdx,3)
                pMin = subPatch % uvRange(patchIdx,1);  pMax = subPatch % uvRange(patchIdx,2)
      case(4);  boundaryVal = subPatch % uvRange(patchIdx,4)
                pMin = subPatch % uvRange(patchIdx,1);  pMax = subPatch % uvRange(patchIdx,2)
    end select

    do i = 1, subPatch % n
      if (i == patchIdx) cycle
      if (subPatch % orig(i) /= P) cycle

      ! Check that Sj shares Si's edge boundary and lies within Si's range
      onEdge = .false.
      select case (edgeIdx)
        case(1)
          onEdge = (abs(subPatch % uvRange(i,2) - boundaryVal) < SURF_TOL) .and. &
                   (subPatch % uvRange(i,3) >= pMin - SURF_TOL) .and. &
                   (subPatch % uvRange(i,4) <= pMax + SURF_TOL)
        case(2)
          onEdge = (abs(subPatch % uvRange(i,1) - boundaryVal) < SURF_TOL) .and. &
                   (subPatch % uvRange(i,3) >= pMin - SURF_TOL) .and. &
                   (subPatch % uvRange(i,4) <= pMax + SURF_TOL)
        case(3)
          onEdge = (abs(subPatch % uvRange(i,4) - boundaryVal) < SURF_TOL) .and. &
                   (subPatch % uvRange(i,1) >= pMin - SURF_TOL) .and. &
                   (subPatch % uvRange(i,2) <= pMax + SURF_TOL)
        case(4)
          onEdge = (abs(subPatch % uvRange(i,3) - boundaryVal) < SURF_TOL) .and. &
                   (subPatch % uvRange(i,1) >= pMin - SURF_TOL) .and. &
                   (subPatch % uvRange(i,2) <= pMax + SURF_TOL)
      end select
      if (.not. onEdge) cycle

      ! Extract the 2 corners of Sj that lie on the shared boundary, and their
      ! true parameter along the edge (p1 for pt1, p2 for pt2).
      ! Si's e=1 (left) is adjacent to Sj's e=2 (right), etc.
      select case (edgeIdx)
        case(1);  pt1 = subPatch % ctrlPts(i,4,1,:);  pt2 = subPatch % ctrlPts(i,4,4,:)
                  p1  = subPatch % uvRange(i,3);       p2  = subPatch % uvRange(i,4)
        case(2);  pt1 = subPatch % ctrlPts(i,1,1,:);  pt2 = subPatch % ctrlPts(i,1,4,:)
                  p1  = subPatch % uvRange(i,3);       p2  = subPatch % uvRange(i,4)
        case(3);  pt1 = subPatch % ctrlPts(i,1,4,:);  pt2 = subPatch % ctrlPts(i,4,4,:)
                  p1  = subPatch % uvRange(i,1);       p2  = subPatch % uvRange(i,2)
        case(4);  pt1 = subPatch % ctrlPts(i,1,1,:);  pt2 = subPatch % ctrlPts(i,4,1,:)
                  p1  = subPatch % uvRange(i,1);       p2  = subPatch % uvRange(i,2)
      end select

      nBufBefore = nBuf
      call insertUnique(buf, nBuf, pt1, SURF_TOL)
      if (nBuf > nBufBefore) paramArr(nBuf) = p1

      nBufBefore = nBuf
      call insertUnique(buf, nBuf, pt2, SURF_TOL)
      if (nBuf > nBufBefore) paramArr(nBuf) = p2
    end do

    if (nBuf < 2) then
      nBuf = 2;  buf(1,:) = ptA;  buf(2,:) = ptB
      return
    end if

    ! Sort by true Bezier parameter along the edge
    call sortByKey(buf, nBuf, paramArr)

  end subroutine collectIntraPatchEdgeVerts

  !!
  !! Collect all unique corner points from sub-patches of adjPatch that lie on
  !! edge adjEdgeIdx of the original patch, sorted along the edge from ptA
  !! toward ptB. Always includes ptA (first) and ptB (last) — buf has at
  !! least 2 entries.
  !!
  !! Sorted by true Bezier parameter (each candidate's own u or v range
  !! boundary, whichever axis varies along this edge), not by projection onto
  !! the straight chord ptA-ptB: chord projection fails whenever the shared
  !! Bezier edge is not monotonic along that chord (e.g. a patch whose control
  !! polygon loops back in one coordinate while bulging out in another),
  !! which can put an interior point before the corner it should follow,
  !! leaving a gap at one segment and a duplicated triangle at another. The
  !! true Bezier parameter is correct regardless of how the physical curve
  !! bends. Direction (whether increasing parameter runs ptA->ptB or
  !! ptB->ptA) is not known a priori -- adjPatch's own parameterisation may
  !! run either way relative to ptA/ptB -- so it is resolved once, after
  !! collection, from the data itself.
  !!
  subroutine collectEdgeVerts(subPatch, adjPatch, adjEdgeIdx, ptA, ptB, buf, nBuf)
    type(subPatches), intent(in)            :: subPatch
    integer(shortInt), intent(in)           :: adjPatch, adjEdgeIdx
    real(defReal), dimension(3), intent(in) :: ptA, ptB
    real(defReal), dimension(:,:), intent(out) :: buf
    integer(shortInt), intent(out)             :: nBuf

    integer(shortInt) :: i, si, nBufBefore
    logical(defBool)  :: onEdge
    real(defReal), dimension(3) :: pt1, pt2
    real(defReal), dimension(size(buf,1)) :: paramArr
    real(defReal) :: p1, p2
    integer(shortInt) :: iMin, iMax

    nBuf = 0

    do i = 1, subPatch % n
      if (subPatch % orig(i) /= adjPatch) cycle

      onEdge = .false.
      select case (adjEdgeIdx)
        case(1);  onEdge = (subPatch % uvRange(i,1) < SURF_TOL)
        case(2);  onEdge = (subPatch % uvRange(i,2) > ONE - SURF_TOL)
        case(3);  onEdge = (subPatch % uvRange(i,3) < SURF_TOL)
        case(4);  onEdge = (subPatch % uvRange(i,4) > ONE - SURF_TOL)
      end select
      if (.not. onEdge) cycle

      select case (adjEdgeIdx)
        case(1);  pt1 = subPatch % ctrlPts(i,1,1,:);  pt2 = subPatch % ctrlPts(i,1,4,:)
                  p1  = subPatch % uvRange(i,3);       p2  = subPatch % uvRange(i,4)
        case(2);  pt1 = subPatch % ctrlPts(i,4,1,:);  pt2 = subPatch % ctrlPts(i,4,4,:)
                  p1  = subPatch % uvRange(i,3);       p2  = subPatch % uvRange(i,4)
        case(3);  pt1 = subPatch % ctrlPts(i,1,1,:);  pt2 = subPatch % ctrlPts(i,4,1,:)
                  p1  = subPatch % uvRange(i,1);       p2  = subPatch % uvRange(i,2)
        case(4);  pt1 = subPatch % ctrlPts(i,1,4,:);  pt2 = subPatch % ctrlPts(i,4,4,:)
                  p1  = subPatch % uvRange(i,1);       p2  = subPatch % uvRange(i,2)
      end select

      nBufBefore = nBuf
      call insertUnique(buf, nBuf, pt1, SURF_TOL)
      if (nBuf > nBufBefore) paramArr(nBuf) = p1

      nBufBefore = nBuf
      call insertUnique(buf, nBuf, pt2, SURF_TOL)
      if (nBuf > nBufBefore) paramArr(nBuf) = p2
    end do

    ! Fallback: if no sub-patches found, just return the two endpoints
    if (nBuf < 2) then
      nBuf = 2;  buf(1,:) = ptA;  buf(2,:) = ptB
      return
    end if

    ! Resolve direction: locate the collected points with the smallest and
    ! largest true parameter (the two true edge endpoints) and check which one
    ! is spatially closer to ptA. If it is the max-parameter one, adjPatch's
    ! parameterisation runs opposite to ptA->ptB, so negate before sorting.
    iMin = 1;  iMax = 1
    do si = 2, nBuf
      if (paramArr(si) < paramArr(iMin)) iMin = si
      if (paramArr(si) > paramArr(iMax)) iMax = si
    end do
    if (norm2(buf(iMax,:) - ptA) < norm2(buf(iMin,:) - ptA)) then
      paramArr(1:nBuf) = -paramArr(1:nBuf)
    end if

    ! Sort by true Bezier parameter along the edge
    call sortByKey(buf, nBuf, paramArr)

  end subroutine collectEdgeVerts

  !!
  !! Merges a neighbouring sub-patch's boundary vertices into this sub-patch's buffer
  !!
  !! Reconciles two adjacent patches' independently-generated T-junction
  !! vertices along a shared edge, when both were selectively subdivided.
  !! collectIntraPatchEdgeVerts only ever looks within one original patch, so
  !! it cannot see a neighbour's independent split points; this merges the
  !! neighbour's boundary vertices (collected separately via collectEdgeVerts)
  !! in, so both sides' triangulations reference an identical vertex set.
  !!
  !! Vertices from extraBuf are accepted only if their projection onto
  !! ptA->ptB falls within [0,1] (extraBuf may come from a neighbour's FULL
  !! boundary edge, which can span more than this sub-patch's own [ptA,ptB]
  !! portion of it, e.g. when the neighbour is subdivided into multiple
  !! pieces along the shared edge and this call only owns one piece). No
  !! perpendicular-distance check is needed: every candidate in extraBuf
  !! already sits exactly on the true shared curve (buildAdjacency verified
  !! the two patches' boundary control points match exactly), so the
  !! projection range is the only thing worth checking. Accepted vertices are
  !! merged into buf with exact-position deduplication, then the whole of buf
  !! is re-sorted by position along ptA->ptB.
  !!
  !! Args:
  !!   buf [inout] -> this sub-patch's edge-vertex buffer; on entry holds
  !!     buf(1:nBuf), on exit holds the merged, re-sorted result
  !!   nBuf [inout] -> number of valid entries in buf
  !!   extraBuf [in] -> the neighbouring sub-patch's edge-vertex buffer to merge in
  !!   nExtra [in] -> number of valid entries in extraBuf
  !!   ptA [in] -> start point of the shared edge, in the direction buf is sorted by
  !!   ptB [in] -> end point of the shared edge
  !!
  !! Result:
  !!   None (buf and nBuf are updated in place).
  !!
  subroutine mergeEdgeVerts(buf, nBuf, extraBuf, nExtra, ptA, ptB)
    real(defReal), dimension(:,:), intent(inout) :: buf
    integer(shortInt), intent(inout)             :: nBuf
    real(defReal), dimension(:,:), intent(in)    :: extraBuf
    integer(shortInt), intent(in)                :: nExtra
    real(defReal), dimension(3), intent(in)      :: ptA, ptB

    integer(shortInt) :: si, sj
    real(defReal), dimension(3) :: edgeDir
    real(defReal), dimension(size(buf,1)) :: dotArr
    real(defReal) :: tParam, lenSq

    edgeDir = ptB - ptA
    lenSq = dot_product(edgeDir, edgeDir)

    do sj = 1, nExtra
      ! Reject candidates outside this sub-patch's [ptA,ptB] span (see doc
      ! comment above).
      if (lenSq > SURF_TOL) then
        tParam = dot_product(extraBuf(sj,:) - ptA, edgeDir) / lenSq
        if (tParam < -SURF_TOL .or. tParam > ONE + SURF_TOL) cycle
      end if

      call insertUnique(buf, nBuf, extraBuf(sj,:), SURF_TOL)
    end do

    ! Re-sort by position along ptA->ptB
    do si = 1, nBuf
      dotArr(si) = dot_product(buf(si,:) - ptA, edgeDir)
    end do
    call sortByKey(buf, nBuf, dotArr)

  end subroutine mergeEdgeVerts

  !!
  !! Checks whether a sub-patch's edge lies on its original patch's outer boundary
  !!
  !! Distinguishes a sub-patch edge that sits on its ORIGINAL patch's true outer
  !! boundary (uMin/uMax/vMin/vMax at 0 or 1) from one purely internal to the
  !! original patch (shared only with another sub-patch of the same original).
  !!
  !! Args:
  !!   uvRange [in] -> UV parameter range table [uMin, uMax, vMin, vMax] for
  !!     every current sub-patch
  !!   idx [in] -> index of the sub-patch to check, into uvRange
  !!   edge [in] -> local edge number to check (1-4; see module header for
  !!     the edge-numbering convention)
  !!
  !! Result:
  !!   .true. if edge `edge` of sub-patch `idx` lies on the [0,1] boundary of
  !!   the original patch's parameter space; .false. if it is an internal
  !!   split introduced by subdivision.
  !!
  pure function edgeOnOuterBoundary(uvRange, idx, edge) result(onBoundary)
    real(defReal), dimension(:,:), intent(in) :: uvRange
    integer(shortInt), intent(in)             :: idx, edge
    logical(defBool)                          :: onBoundary

    onBoundary = .false.
    select case (edge)
      case(1);  onBoundary = uvRange(idx,1) < SURF_TOL
      case(2);  onBoundary = uvRange(idx,2) > ONE - SURF_TOL
      case(3);  onBoundary = uvRange(idx,3) < SURF_TOL
      case(4);  onBoundary = uvRange(idx,4) > ONE - SURF_TOL
    end select

  end function edgeOnOuterBoundary

  !!
  !! Build a polygon from two edge-vertex buffers sharing a corner (bufA in
  !! full, then bufB's interior points bufB(2:nB)), and fan-triangulate it
  !! from that shared corner (poly(1)) into tris, skipping degenerate
  !! triangles.
  !!
  !! Args:
  !!   bufA [in], nA [in] -> first edge buffer, included in full
  !!   bufB [in], nB [in] -> second edge buffer, only points 2:nB appended
  !!   tris [inout] -> triangle mesh being built
  !!   nTri [inout] -> number of valid triangles in tris
  !!
  !! Errors:
  !!   fatalError if the combined polygon would exceed MAX_POLY.
  !!
  subroutine fanTriangulate(bufA, nA, bufB, nB, tris, nTri)
    real(defReal), dimension(:,:), intent(in)      :: bufA, bufB
    integer(shortInt), intent(in)                  :: nA, nB
    real(defReal), dimension(:,:,:), intent(inout) :: tris
    integer(shortInt), intent(inout)               :: nTri
    real(defReal), dimension(MAX_POLY, 3) :: poly
    integer(shortInt) :: si, nPoly
    character(100), parameter :: Here = 'fanTriangulate (bezierVolume_class.f90)'

    nPoly = nA
    poly(1:nA, :) = bufA(1:nA, :)
    do si = 2, nB
      nPoly = nPoly + 1
      if (nPoly > MAX_POLY) call fatalError(Here, 'MAX_POLY exceeded ('//numToChar(MAX_POLY)//')')
      poly(nPoly, :) = bufB(si, :)
    end do
    do si = 2, nPoly - 1
      if (.not. triDegenerate(poly(1,:), poly(si,:), poly(si+1,:))) then
        nTri = nTri + 1
        tris(nTri, 1, :) = poly(1, :)
        tris(nTri, 2, :) = poly(si, :)
        tris(nTri, 3, :) = poly(si+1, :)
      end if
    end do

  end subroutine fanTriangulate

  !!
  !! True if triangle (v0,v1,v2) has two coincident vertices (zero area).
  !!
  !! A fan-cap patch whose whole edge1 row is collapsed to a single point
  !! (e.g. a polar apex) has C00 == C01 exactly. The unsubdivided-patch T1 fan
  !! (poly = e1buf || e4buf[2:], fan from C00) then emits a triangle
  !! (C00, C00, C11) with two identical vertices -- zero area, contributes
  !! nothing to any real ray cast, but its two "edges" that aren't truly zero
  !! length duplicate T2's genuine edges, showing up as spurious duplicate
  !! edges in a watertightness check. Filtering these out at emission is a
  !! strict improvement: a degenerate triangle can never be a needed hit.
  !!
  pure function triDegenerate(v0, v1, v2) result(deg)
    real(defReal), dimension(3), intent(in) :: v0, v1, v2
    logical(defBool) :: deg
    deg = (norm2(v0-v1) < SURF_TOL) .or. (norm2(v1-v2) < SURF_TOL) .or. (norm2(v0-v2) < SURF_TOL)
  end function triDegenerate

  !!
  !! Ray cast against the provided triangle list to determine inside/outside.
  !! Counts forward intersections; odd = inside, even = outside.
  !!
  function rayCast(self, tris, nTri, r) result(inside)
    class(bezierVolume), intent(in)              :: self
    real(defReal), dimension(:,:,:), intent(in)  :: tris     ! (nTri, 3, 3)
    integer(shortInt), intent(in)                :: nTri
    real(defReal), dimension(3), intent(in)      :: r
    logical(defBool)                             :: inside
    real(defReal), dimension(3) :: rayDir, v0, v1, v2, rLocal
    integer(shortInt)           :: count, i
    logical(defBool)            :: hit, nearZero, anyNearZero

    rayDir = (/ ONE, ONE / 3.0_defReal, ONE / 7.0_defReal /)
    rayDir = rayDir / norm2(rayDir)

    count = 0
    anyNearZero = .false.

    do i = 1, nTri
      v0 = tris(i, 1, :);  v1 = tris(i, 2, :);  v2 = tris(i, 3, :)
      call self % rayTriangle(r, rayDir, v0, v1, v2, hit, nearZero)
      if (hit)     count = count + 1
      if (nearZero) anyNearZero = .true.
    end do

    ! If degenerate and no hits, nudge the origin and retry
    if (count == 0 .and. anyNearZero) then
      rLocal = r + (/ NUDGE, NUDGE, NUDGE /)
      count = 0
      do i = 1, nTri
        v0 = tris(i, 1, :);  v1 = tris(i, 2, :);  v2 = tris(i, 3, :)
        call self % rayTriangle(rLocal, rayDir, v0, v1, v2, hit, nearZero)
        if (hit) count = count + 1
      end do
    end if

    inside = mod(count, 2) == 1

  end function rayCast

  !!
  !! Watertight ray-triangle intersection (Woop, Benthin & Wald 2013).
  !! Edge functions depend only on each edge's two vertices, so shared edges
  !! between adjacent triangles produce bitwise-identical results (no gaps).
  !!
  !! Args:
  !!   self [in] -> the bezierVolume (unused; part of the type-bound interface)
  !!   r [in] -> ray origin
  !!   rayDir [in] -> ray direction (need not be normalised)
  !!   v0 [in], v1 [in], v2 [in] -> the triangle's 3 vertices
  !!   hit [out] -> .true. if the ray hits the triangle at a positive distance
  !!   nearZero [out] -> .true. if the hit distance is close to zero relative
  !!     to the triangle's own scale (rayCast's retry-with-nudged-origin
  !!     trigger for query points sitting almost exactly on a triangle)
  !!
  subroutine rayTriangle(self, r, rayDir, v0, v1, v2, hit, nearZero)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, rayDir, v0, v1, v2
    logical(defBool), intent(out)           :: hit, nearZero
    real(defReal), dimension(3) :: A, B, C
    real(defReal) :: Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz
    real(defReal) :: e0, e1, e2, det, t
    real(defReal) :: abs1, abs2, abs3, Sx, Sy, Sz
    integer(shortInt) :: kz, kx, ky
    real(defReal), parameter :: NEAR_ZERO_TOL = 1.0E-4_defReal

    hit      = .false.
    nearZero = .false.

    abs1 = abs(rayDir(1))
    abs2 = abs(rayDir(2))
    abs3 = abs(rayDir(3))

    if (abs1 > abs2 .and. abs1 > abs3) then
      kz = 1;  kx = 2;  ky = 3
    else if (abs2 > abs3) then
      kz = 2;  kx = 3;  ky = 1
    else
      kz = 3;  kx = 1;  ky = 2
    end if

    Sx = rayDir(kx) / rayDir(kz)
    Sy = rayDir(ky) / rayDir(kz)
    Sz = ONE        / rayDir(kz)

    A = v0 - r;  B = v1 - r;  C = v2 - r

    Ax = A(kx) - Sx * A(kz);  Ay = A(ky) - Sy * A(kz)
    Bx = B(kx) - Sx * B(kz);  By = B(ky) - Sy * B(kz)
    Cx = C(kx) - Sx * C(kz);  Cy = C(ky) - Sy * C(kz)

    e0 = Bx * Cy - By * Cx
    e1 = Cx * Ay - Cy * Ax
    e2 = Ax * By - Ay * Bx

    if (e0 < -SURF_TOL .or. e1 < -SURF_TOL .or. e2 < -SURF_TOL) then
      if (e0 > SURF_TOL .or. e1 > SURF_TOL .or. e2 > SURF_TOL) return
    end if

    det = e0 + e1 + e2
    if (abs(det) < SURF_TOL) then
      nearZero = .true.
      return
    end if

    Az = Sz * A(kz);  Bz = Sz * B(kz);  Cz = Sz * C(kz)
    t  = e0 * Az + e1 * Bz + e2 * Cz

    if (abs(t) < NEAR_ZERO_TOL * abs(det)) nearZero = .true.

    if (det > ZERO) then
      if (t < SURF_TOL * det) return
    else
      if (t > SURF_TOL * det) return
    end if

    hit = .true.

  end subroutine rayTriangle

  ! ---------------------------------------------------------------------------
  ! GJK convex hull containment test
  ! ---------------------------------------------------------------------------

  !!
  !! GJK test: is point r inside the convex hull of pts(1:nPts,:)?
  !!
  !! Builds a simplex (up to 4 points) by repeatedly finding the support
  !! point furthest from r in the current search direction d (gjkSupport)
  !! and folding it into the simplex (gjkDoSimplex), which either reports the
  !! origin (i.e. r) enclosed or narrows d toward it. Terminates early with
  !! inside = .false. as soon as a support point fails to pass the origin's
  !! side (proof that no further progress toward r is possible); if MAX_ITER
  !! is reached without resolving either way, conservatively reports inside.
  !!
  !! Args:
  !!   pts [in], nPts [in] -> candidate points whose convex hull is tested
  !!   r [in] -> the point being tested
  !!
  !! Result:
  !!   .true. if r is inside (or on the boundary of) the convex hull of
  !!   pts(1:nPts,:), .false. otherwise
  !!
  pure function pointInConvexHull(pts, nPts, r) result(inside)
    real(defReal), dimension(:,:), intent(in) :: pts
    integer(shortInt), intent(in)             :: nPts
    real(defReal), dimension(3), intent(in)   :: r
    logical(defBool)                          :: inside
    real(defReal), dimension(4, 3) :: S
    integer(shortInt)              :: nS, i, iter
    real(defReal), dimension(3)    :: d, sup
    logical(defBool)               :: ok
    integer(shortInt), parameter   :: MAX_ITER = 64

    inside = .false.
    if (nPts < 1) return

    ! Initial direction: from r toward centroid of the point set
    d = ZERO
    do i = 1, nPts
      d = d + pts(i,:)
    end do
    d = d / real(nPts, defReal) - r

    if (dot_product(d, d) < SURF_TOL) then
      inside = .true.
      return
    end if

    call gjkSupport(pts, nPts, r, d, sup, ok)
    if (.not. ok) return

    nS = 1
    S(1,:) = sup
    d = -sup

    do iter = 1, MAX_ITER
      call gjkSupport(pts, nPts, r, d, sup, ok)
      if (.not. ok) return

      nS = nS + 1
      S(nS,:) = sup

      call gjkDoSimplex(S, nS, d, inside)
      if (inside) return
    end do

    inside = .true.  ! failed to converge: conservatively inside

  end function pointInConvexHull

  !!
  !! GJK support function: the point of pts(1:nPts,:) (relative to r) that is
  !! furthest in direction d.
  !!
  !! Args:
  !!   pts [in], nPts [in] -> candidate points
  !!   r [in] -> the point being tested, i.e. the GJK search origin
  !!   d [in] -> search direction
  !!   sup [out] -> the support point, relative to r
  !!   ok [out] -> .false. if sup does not pass the origin's side of the
  !!     plane through it perpendicular to d -- i.e. direction d cannot
  !!     reach the origin from here, so the hull cannot contain r.
  !!     .true. otherwise.
  !!
  pure subroutine gjkSupport(pts, nPts, r, d, sup, ok)
    real(defReal), dimension(:,:), intent(in) :: pts
    integer(shortInt), intent(in)             :: nPts
    real(defReal), dimension(3), intent(in)   :: r, d
    real(defReal), dimension(3), intent(out)  :: sup
    logical(defBool), intent(out)             :: ok
    integer(shortInt) :: i, bestIdx
    real(defReal) :: maxDot, dp

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
    ok = (dot_product(sup, d) >= ZERO)

  end subroutine gjkSupport

  !!
  !! GJK simplex-evolution step: given the current simplex S (nS points) and
  !! search direction d, either detect that the simplex encloses the origin
  !! (inside = .true.) or discard the point(s) not facing the origin and
  !! update d to the next search direction. Dispatches to gjkTetrahedron,
  !! gjkTriangle or gjkLine depending on the simplex size.
  !!
  !! Args:
  !!   S [inout] -> simplex points, S(1:nS,3); reduced/reordered in place
  !!   nS [inout] -> number of valid points in S; reduced in place
  !!   d [inout] -> search direction; updated in place unless inside
  !!   inside [out] -> true once the simplex is found to enclose the origin
  !!
  pure subroutine gjkDoSimplex(S, nS, d, inside)
    real(defReal), dimension(4,3), intent(inout) :: S
    integer(shortInt), intent(inout)             :: nS
    real(defReal), dimension(3), intent(inout)   :: d
    logical(defBool), intent(out)                :: inside

    inside = .false.

    select case (nS)
      case (4)
        call gjkTetrahedron(S, nS, d, inside)
      case (3)
        call gjkTriangle(S, nS, d, inside)
      case (2)
        call gjkLine(S, nS, d)
    end select

  end subroutine gjkDoSimplex

  !!
  !! Reduce a 4-point tetrahedron simplex S(1:4,:) (origin vertex a = S(4,:)).
  !! If the tetrahedron encloses the origin, reports inside = .true.
  !! Otherwise, discards the face facing away from the origin and continues
  !! via gjkTriangle on the remaining 3-point face (explicit fall-through,
  !! rather than the implicit fall-through a shared nS would otherwise give).
  !!
  pure subroutine gjkTetrahedron(S, nS, d, inside)
    real(defReal), dimension(4,3), intent(inout) :: S
    integer(shortInt), intent(inout)             :: nS
    real(defReal), dimension(3), intent(inout)   :: d
    logical(defBool), intent(out)                :: inside
    real(defReal), dimension(3) :: a, b, c, dNext
    real(defReal), dimension(3) :: ab, ac, ad, ao
    real(defReal), dimension(3) :: abc, acd, adb

    a  = S(4,:);  b = S(3,:);  c = S(2,:);  dNext = S(1,:)
    ab = b - a;   ac = c - a;  ad = dNext - a;  ao = -a

    abc = cross3(ab, ac);  if (dot_product(abc, ad) > ZERO) abc = -abc
    acd = cross3(ac, ad);  if (dot_product(acd, ab) > ZERO) acd = -acd
    adb = cross3(ad, ab);  if (dot_product(adb, ac) > ZERO) adb = -adb

    if (dot_product(abc, ao) > ZERO) then
      if (dot_product(cross3(ab, ac), ad) < ZERO) then
        S(1,:) = c;  S(2,:) = b;  S(3,:) = a
      else
        S(1,:) = b;  S(2,:) = c;  S(3,:) = a
      end if
    else if (dot_product(acd, ao) > ZERO) then
      if (dot_product(cross3(ac, ad), ab) < ZERO) then
        S(1,:) = dNext;  S(2,:) = c;  S(3,:) = a
      else
        S(1,:) = c;  S(2,:) = dNext;  S(3,:) = a
      end if
    else if (dot_product(adb, ao) > ZERO) then
      if (dot_product(cross3(ad, ab), ac) < ZERO) then
        S(1,:) = b;  S(2,:) = dNext;  S(3,:) = a
      else
        S(1,:) = dNext;  S(2,:) = b;  S(3,:) = a
      end if
    else
      inside = .true.
      return
    end if

    nS = 3
    call gjkTriangle(S, nS, d, inside)

  end subroutine gjkTetrahedron

  !!
  !! Reduce a 3-point triangle simplex S(1:3,:) (origin vertex a = S(3,:)).
  !! Never finds the origin enclosed (a 2D triangle cannot enclose a 3D
  !! origin) -- always either discards a vertex/edge, reducing to a 2-point
  !! or 1-point simplex, or keeps the triangle and points d at whichever side
  !! of it faces the origin.
  !!
  pure subroutine gjkTriangle(S, nS, d, inside)
    real(defReal), dimension(4,3), intent(inout) :: S
    integer(shortInt), intent(inout)             :: nS
    real(defReal), dimension(3), intent(inout)   :: d
    logical(defBool), intent(out)                :: inside
    real(defReal), dimension(3) :: a, b, c
    real(defReal), dimension(3) :: ab, ac, ao
    real(defReal), dimension(3) :: abc, abPerp, acPerp
    logical(defBool)            :: checkAB

    inside = .false.
    a  = S(3,:);  b = S(2,:);  c = S(1,:)
    ab = b - a;   ac = c - a;  ao = -a

    abc    = cross3(ab, ac)
    acPerp = cross3(abc, ac)
    abPerp = cross3(ab, abc)
    checkAB = .false.

    if (dot_product(acPerp, ao) > ZERO) then
      if (dot_product(ac, ao) > ZERO) then
        S(1,:) = c;  S(2,:) = a;  nS = 2
        d = tripleCross(ac, ao, ac)
        return
      else
        checkAB = .true.
      end if
    end if

    if (checkAB .or. dot_product(abPerp, ao) > ZERO) then
      if (dot_product(ab, ao) > ZERO) then
        S(1,:) = b;  S(2,:) = a;  nS = 2
        d = tripleCross(ab, ao, ab)
      else
        S(1,:) = a;  nS = 1;  d = ao
      end if
      return
    end if

    if (dot_product(abc, ao) > ZERO) then
      d = abc
    else
      S(1,:) = b;  S(2,:) = c;  S(3,:) = a
      d = -abc
    end if

  end subroutine gjkTriangle

  !!
  !! Reduce a 2-point line-segment simplex S(1:2,:) (origin vertex a = S(2,:)).
  !! Either discards a vertex, reducing to a 1-point simplex pointing straight
  !! at the origin, or keeps the line and points d perpendicular to it,
  !! toward the origin (falling back to an arbitrary perpendicular if the
  !! origin lies exactly on the line).
  !!
  pure subroutine gjkLine(S, nS, d)
    real(defReal), dimension(4,3), intent(inout) :: S
    integer(shortInt), intent(inout)             :: nS
    real(defReal), dimension(3), intent(inout)   :: d
    real(defReal), dimension(3) :: a, b, ab, ao, tcross

    a  = S(2,:);  b = S(1,:)
    ab = b - a;   ao = -a

    if (dot_product(ab, ao) > ZERO) then
      tcross = tripleCross(ab, ao, ab)
      if (dot_product(tcross, tcross) < SURF_TOL) then
        if (abs(ab(1)) < abs(ab(2))) then
          d(1) = ZERO;  d(2) = -ab(3);  d(3) = ab(2)
        else
          d(1) = ab(3);  d(2) = ZERO;  d(3) = -ab(1)
        end if
      else
        d = tcross
      end if
    else
      S(1,:) = a;  nS = 1;  d = ao
    end if

  end subroutine gjkLine

  !!
  !! Standard right-handed cross product a x b.
  !!
  !! NOTE: this is NOT the same sign convention as genericProcedures::crossProduct,
  !! which returns b x a (the negated result). The two are not interchangeable --
  !! GJK's simplex-winding logic here depends on the standard a x b sign.
  !!
  pure function cross3(a, b) result(c)
    real(defReal), dimension(3), intent(in) :: a, b
    real(defReal), dimension(3)             :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  !!
  !! Triple cross product (a x b) x c.
  !!
  pure function tripleCross(a, b, c) result(r)
    real(defReal), dimension(3), intent(in) :: a, b, c
    real(defReal), dimension(3)             :: r
    r = cross3(cross3(a, b), c)
  end function tripleCross

end module bezierVolume_class
