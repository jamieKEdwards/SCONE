module bezierVolume_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill
  implicit none
  private

  ! Counter for particles that hit the subdivision limit
  integer(shortInt), public :: bezierVolumeExceededCount = 0

  ! Safety limit on subdivision depth per halfspace call
  integer(shortInt), parameter :: MAX_SUBDIVISIONS = 10

  ! GJK tolerance
  real(defReal), parameter :: GJK_EPS = 1.0E-12_defReal

  !!
  !! 3D Bezier volume surface defined by a watertight set of bicubic Bezier patches.
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
    procedure :: halfspace
    procedure :: inPatchAABB
    procedure :: rayCast
    procedure :: rayCastDebug
    procedure :: rayTriangle
    procedure :: buildAdjacency
  end type bezierVolume


contains

  ! ---------------------------------------------------------------------------
  ! Surface interface routines
  ! ---------------------------------------------------------------------------

  pure function myType(self) result(str)
    class(bezierVolume), intent(in) :: self
    character(:), allocatable       :: str
    str = 'bezierVolume'
  end function myType

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
  !! Build the patch adjacency map self%adj(numPatches, 4).
  !! Two patches share an edge if all 4 boundary control points match
  !! (forward or reverse order) to within MATCH_TOL.
  !!
  subroutine buildAdjacency(self)
    class(bezierVolume), intent(inout) :: self
    integer(shortInt) :: i, j, ei, ej, k
    real(defReal), dimension(4,3) :: edgeI, edgeJ
    logical(defBool) :: fwd, rev
    real(defReal), parameter :: MATCH_TOL = 1.0E-10_defReal

    allocate(self % adj(self % numPatches, 4))
    allocate(self % adjEdge(self % numPatches, 4))
    self % adj     = -1
    self % adjEdge = -1

    do i = 1, self % numPatches
      do ei = 1, 4
        edgeI = patchEdge(self % ctrlPts(i,:,:,:), ei)
        do j = i + 1, self % numPatches
          do ej = 1, 4
            edgeJ = patchEdge(self % ctrlPts(j,:,:,:), ej)
            fwd = .true.
            rev = .true.
            do k = 1, 4
              if (any(abs(edgeI(k,:) - edgeJ(k,:))   > MATCH_TOL)) fwd = .false.
              if (any(abs(edgeI(k,:) - edgeJ(5-k,:)) > MATCH_TOL)) rev = .false.
            end do
            if (fwd .or. rev) then
              self % adj(i, ei)     = j;  self % adjEdge(i, ei) = ej
              self % adj(j, ej)     = i;  self % adjEdge(j, ej) = ei
            end if
          end do
        end do
      end do
    end do

  end subroutine buildAdjacency

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

  pure function evaluate(self, r) result(c)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    c = ZERO
  end function evaluate

  pure function distance(self, r, u) result(d)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: d
    d = INF
  end function distance

  pure function going(self, r, u) result(halfspace)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: halfspace
    halfspace = .false.
  end function going

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
  !! Determine halfspace for a particle position.
  !!
  !! Result: .true. = outside (positive halfspace)
  !!
  function halfspace(self, r, u) result(hs)
    class(bezierVolume), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    logical(defBool)                        :: hs

    ! Selective subdivision working arrays
    real(defReal), dimension(:,:,:,:), allocatable :: curPts
    real(defReal), dimension(:,:,:),   allocatable :: curWeights, newWeights
    integer(shortInt), dimension(:), allocatable   :: origPatch, newOrig
    logical(defBool), dimension(:), allocatable    :: needsSubdiv
    real(defReal), dimension(:,:,:,:), allocatable :: newPts
    ! UV parameter ranges for each sub-patch: [uMin, uMax, vMin, vMax]
    real(defReal), dimension(:,:), allocatable :: uvRange, newUV

    ! Subdivision tracking
    logical(defBool), dimension(self % numPatches) :: wasSubdivided
    ! Original patches that have passed the STRICT hull test at least once this
    ! call. Used to scope the AABB-only fallback (see Step 3) to patches actually
    ! adjacent to the region of interest, instead of scanning every patch in the
    ! model by loose bounding-box coincidence.
    logical(defBool), dimension(self % numPatches) :: everHullSelected
    ! Candidate original-patch index for the AABB-fallback eligibility check
    integer(shortInt) :: candidatePatch

    ! Triangle mesh for ray cast
    real(defReal), dimension(:,:,:), allocatable :: tris
    integer(shortInt) :: nTri, maxTri

    ! Working variables
    real(defReal), dimension(16, 3) :: patchPts
    real(defReal), dimension(4,4,3) :: Q00, Q01, Q10, Q11
    real(defReal), dimension(4,4)   :: Q00w, Q01w, Q10w, Q11w
    integer(shortInt) :: subdivision, i, j, k, m, nCur, newN, nNeed, P
    integer(shortInt) :: eP_idx, P_adj

    ! Corner points for mesh building
    real(defReal), dimension(3) :: C00, C01, C10, C11

    ! T-junction vertex collection: per-edge buffers + sort keys
    ! MAX_TJ_BUF handles up to 2^8 = 256 sub-patches along one edge
    integer(shortInt), parameter :: MAX_TJ_BUF = 260
    real(defReal), dimension(MAX_TJ_BUF, 3) :: e1buf, e2buf, e3buf, e4buf
    real(defReal), dimension(MAX_TJ_BUF)    :: e1dot, e2dot, e3dot, e4dot
    integer(shortInt) :: n1, n2, n3, n4
    ! Cross-patch merge buffer: when a subdivided patch's edge is on the ORIGINAL
    ! patch's outer boundary and the adjacent original patch was ALSO subdivided,
    ! neither side's independent quadtree refinement is synchronised with the
    ! other along their shared edge. collectIntraPatchEdgeVerts only ever looks
    ! within the same original patch, so it cannot see the neighbour's split
    ! points. mergeEdgeVerts pulls the neighbour's boundary vertices in via
    ! collectEdgeVerts and merges them into the same buffer used for the fan,
    ! so both sides' triangulations reference an identical vertex set.
    real(defReal), dimension(MAX_TJ_BUF, 3) :: xBuf
    integer(shortInt) :: nX
    logical(defBool) :: onEdge, isDup
    real(defReal), dimension(3) :: pt1, pt2, tmpV3
    real(defReal) :: tmpDot, u0, u1, v0, v1, um, vm
    integer(shortInt) :: si
    real(defReal), parameter :: TJ_TOL = 1.0E-10_defReal

    ! Polygon buffer for fan triangulation of each triangle face
    integer(shortInt), parameter :: MAX_POLY = 530   ! 2 * MAX_TJ_BUF + endpoints
    real(defReal), dimension(MAX_POLY, 3) :: poly
    integer(shortInt) :: nPoly

    ! --- Step 1: AABB rejection ---
    if (r(1) < self % aabb(1) - SURF_TOL .or. r(1) > self % aabb(4) + SURF_TOL .or. &
        r(2) < self % aabb(2) - SURF_TOL .or. r(2) > self % aabb(5) + SURF_TOL .or. &
        r(3) < self % aabb(3) - SURF_TOL .or. r(3) > self % aabb(6) + SURF_TOL) then
      hs = .true.
      go to 999
    end if

    ! --- Step 2: Global convex hull rejection ---
    if (.not. pointInConvexHull(self % allPtsFlat, self % nAllPts, r)) then
      hs = .true.
      go to 999
    end if

    ! --- Step 3: Selective subdivision ---
    nCur = self % numPatches
    allocate(curPts(nCur, 4, 4, 3))
    allocate(curWeights(nCur, 4, 4))
    allocate(origPatch(nCur))
    allocate(uvRange(nCur, 4))
    curPts     = self % ctrlPts
    curWeights = self % weights
    do i = 1, nCur
      origPatch(i) = i
      uvRange(i, :) = (/ ZERO, ONE, ZERO, ONE /)
    end do
    wasSubdivided = .false.
    everHullSelected = .false.

    do subdivision = 1, MAX_SUBDIVISIONS

      ! Identify which current patches need subdivision (AABB + GJK on 16 control points)
      allocate(needsSubdiv(nCur))
      needsSubdiv = .false.
      nNeed = 0

      do i = 1, nCur
        if (.not. self % inPatchAABB(curPts(i,:,:,:), r)) cycle
        m = 0
        do j = 1, 4
          do k = 1, 4
            m = m + 1
            patchPts(m,:) = curPts(i, j, k, :)
          end do
        end do
        if (pointInConvexHull(patchPts, 16, r)) then
          needsSubdiv(i) = .true.
          nNeed = nNeed + 1
          everHullSelected(origPatch(i)) = .true.
        end if
      end do

      if (nNeed == 0) then
        ! GJK found no patch: the query point lies on a patch boundary where no
        ! single patch's convex hull strictly contains it (typically floating-point
        ! precision on an already-tight, deeply-subdivided hull). Fall back to
        ! AABB-only selection so boundary points still trigger subdivision.
        !
        ! Scope: if any original patch has EVER passed the strict hull test this
        ! call, restrict the fallback to patches that are that patch (still-live
        ! sub-patches of it) or genuinely adjacent to it -- not every patch in the
        ! model. An unscoped scan pulls in totally unrelated patches purely by
        ! loose bounding-box coincidence (confirmed 2026-08-01 on Gumbo: patch 15
        ! was pulled into subdivision this way despite its own hull never once
        ! containing the query point, then only ever refined one level deep while
        ! its genuine neighbour kept refining much further -- a mismatch that
        ! opened a real gap between them; see bezierVolume_status.md).
        ! On the first round with no prior hull match at all, there is no adjacency
        ! information yet to scope by, so fall back to the original unrestricted
        ! scan (matches pre-fix behaviour for that narrow case).
        do i = 1, nCur
          if (.not. self % inPatchAABB(curPts(i,:,:,:), r)) cycle
          candidatePatch = origPatch(i)
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
      newN = nCur + 3 * nNeed
      allocate(newPts(newN, 4, 4, 3))
      allocate(newWeights(newN, 4, 4))
      allocate(newOrig(newN))
      allocate(newUV(newN, 4))
      j = 0
      do i = 1, nCur
        if (needsSubdiv(i)) then
          wasSubdivided(origPatch(i)) = .true.
          call subdividePatch(curPts(i,:,:,:), curWeights(i,:,:), &
                              Q00, Q01, Q10, Q11, Q00w, Q01w, Q10w, Q11w)
          u0 = uvRange(i,1);  u1 = uvRange(i,2)
          v0 = uvRange(i,3);  v1 = uvRange(i,4)
          um = HALF*(u0+u1);  vm = HALF*(v0+v1)
          ! Q00: u∈[u0,um], v∈[v0,vm]
          j = j + 1;  newPts(j,:,:,:) = Q00;  newWeights(j,:,:) = Q00w
          newOrig(j) = origPatch(i);  newUV(j,:) = (/ u0, um, v0, vm /)
          ! Q01: u∈[u0,um], v∈[vm,v1]
          j = j + 1;  newPts(j,:,:,:) = Q01;  newWeights(j,:,:) = Q01w
          newOrig(j) = origPatch(i);  newUV(j,:) = (/ u0, um, vm, v1 /)
          ! Q10: u∈[um,u1], v∈[v0,vm]
          j = j + 1;  newPts(j,:,:,:) = Q10;  newWeights(j,:,:) = Q10w
          newOrig(j) = origPatch(i);  newUV(j,:) = (/ um, u1, v0, vm /)
          ! Q11: u∈[um,u1], v∈[vm,v1]
          j = j + 1;  newPts(j,:,:,:) = Q11;  newWeights(j,:,:) = Q11w
          newOrig(j) = origPatch(i);  newUV(j,:) = (/ um, u1, vm, v1 /)
        else
          j = j + 1
          newPts(j,:,:,:)  = curPts(i,:,:,:)
          newWeights(j,:,:) = curWeights(i,:,:)
          newOrig(j)        = origPatch(i)
          newUV(j,:)        = uvRange(i,:)
        end if
      end do

      deallocate(needsSubdiv)
      call move_alloc(newPts,     curPts)
      call move_alloc(newWeights, curWeights)
      call move_alloc(newOrig,    origPatch)
      call move_alloc(newUV,      uvRange)
      nCur = newN

    end do

    if (subdivision > MAX_SUBDIVISIONS) then
      !$omp atomic
      bezierVolumeExceededCount = bezierVolumeExceededCount + 1
    end if

    ! --- Step 4: Build watertight triangle mesh with T-junction patching ---
    !
    ! For each original patch P:
    !   If P was subdivided: add 2 corner triangles per sub-patch (watertight by
    !     construction since all sub-patches of P share vertices at their shared edges).
    !   If P was NOT subdivided: add the standard 2 corner triangles, but for any
    !     edge where the adjacent patch WAS subdivided, insert the De Casteljau midpoint
    !     of that boundary edge (T-junction vertex) to split the affected triangle.
    !
    ! Triangle orientation (winding unused — we just count hits):
    !   T1 = (C00, C01, C11)   covers the u=0 and v=1 boundaries
    !   T2 = (C00, C11, C10)   covers the u=1 and v=0 boundaries
    !
    ! Per-patch triangle splits when adjacent patch is subdivided:
    !   Edge 1 (u=0, C00-C01) in T1: M1 on C00-C01  -> split T1 at M1
    !   Edge 4 (v=1, C01-C11) in T1: M4 on C01-C11  -> split T1 at M4
    !   Edge 2 (u=1, C10-C11) in T2: M2 on C10-C11  -> split T2 at M2
    !   Edge 3 (v=0, C00-C10) in T2: M3 on C00-C10  -> split T2 at M3

    maxTri = MAX_TJ_BUF * nCur + MAX_TJ_BUF * self % numPatches
    allocate(tris(maxTri, 3, 3))
    nTri = 0

    do P = 1, self % numPatches

      if (wasSubdivided(P)) then
        ! Polygon-fan triangulation for each sub-patch of P, with within-patch
        ! T-junction detection. Selective subdivision may leave adjacent sub-patches
        ! of the same original patch at different depths; this closes those gaps.
        do i = 1, nCur
          if (origPatch(i) /= P) cycle
          C00 = curPts(i, 1, 1, :)
          C01 = curPts(i, 1, 4, :)
          C10 = curPts(i, 4, 1, :)
          C11 = curPts(i, 4, 4, :)

          call collectIntraPatchEdgeVerts(curPts, nCur, origPatch, uvRange, i, 1, &
                                          C00, C01, e1buf, n1)
          call collectIntraPatchEdgeVerts(curPts, nCur, origPatch, uvRange, i, 2, &
                                          C10, C11, e2buf, n2)
          call collectIntraPatchEdgeVerts(curPts, nCur, origPatch, uvRange, i, 3, &
                                          C00, C10, e3buf, n3)
          call collectIntraPatchEdgeVerts(curPts, nCur, origPatch, uvRange, i, 4, &
                                          C01, C11, e4buf, n4)

          ! Reconcile against an ALSO-subdivided neighbour along P's true outer
          ! boundary edges (see xBuf declaration comment above). Skipped for edges
          ! internal to P (shared with another sub-patch of the same original
          ! patch) -- those are already handled by collectIntraPatchEdgeVerts above.
          if (edgeOnOuterBoundary(uvRange, i, 1)) then
            P_adj = self % adj(P, 1)
            if (P_adj > 0) then
              if (wasSubdivided(P_adj)) then
                eP_idx = self % adjEdge(P, 1)
                call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                      C00, C01, xBuf, nX)
                call mergeEdgeVerts(e1buf, n1, xBuf, nX, C00, C01)
              end if
            end if
          end if
          if (edgeOnOuterBoundary(uvRange, i, 2)) then
            P_adj = self % adj(P, 2)
            if (P_adj > 0) then
              if (wasSubdivided(P_adj)) then
                eP_idx = self % adjEdge(P, 2)
                call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                      C10, C11, xBuf, nX)
                call mergeEdgeVerts(e2buf, n2, xBuf, nX, C10, C11)
              end if
            end if
          end if
          if (edgeOnOuterBoundary(uvRange, i, 3)) then
            P_adj = self % adj(P, 3)
            if (P_adj > 0) then
              if (wasSubdivided(P_adj)) then
                eP_idx = self % adjEdge(P, 3)
                call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                      C00, C10, xBuf, nX)
                call mergeEdgeVerts(e3buf, n3, xBuf, nX, C00, C10)
              end if
            end if
          end if
          if (edgeOnOuterBoundary(uvRange, i, 4)) then
            P_adj = self % adj(P, 4)
            if (P_adj > 0) then
              if (wasSubdivided(P_adj)) then
                eP_idx = self % adjEdge(P, 4)
                call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                      C01, C11, xBuf, nX)
                call mergeEdgeVerts(e4buf, n4, xBuf, nX, C01, C11)
              end if
            end if
          end if

          ! T1: fan over polygon C00 → [e1 interior] → C01 → [e4 interior] → C11
          nPoly = n1
          poly(1:n1, :) = e1buf(1:n1, :)
          do si = 2, n4
            nPoly = nPoly + 1
            poly(nPoly, :) = e4buf(si, :)
          end do
          do si = 2, nPoly - 1
            nTri = nTri + 1
            tris(nTri, 1, :) = poly(1, :)
            tris(nTri, 2, :) = poly(si, :)
            tris(nTri, 3, :) = poly(si+1, :)
          end do

          ! T2: fan over polygon C00 → [e3 interior] → C10 → [e2 interior] → C11
          nPoly = n3
          poly(1:n3, :) = e3buf(1:n3, :)
          do si = 2, n2
            nPoly = nPoly + 1
            poly(nPoly, :) = e2buf(si, :)
          end do
          do si = 2, nPoly - 1
            nTri = nTri + 1
            tris(nTri, 1, :) = poly(1, :)
            tris(nTri, 2, :) = poly(si, :)
            tris(nTri, 3, :) = poly(si+1, :)
          end do
        end do

      else
        ! Unsubdivided patch Q: polygon-fan triangulation with full T-junction matching.
        !
        ! For each boundary edge of Q adjacent to a subdivided patch, collect ALL
        ! sub-patch corners lying on that shared edge, sort them along the edge, and
        ! build a fan from the opposite vertex.  This closes multi-level T-junction
        ! gaps (not just the single t=0.5 midpoint as before).
        !
        ! T1 polygon: C00 → [e1 interior] → C01 → [e4 interior] → C11  (fan from C00)
        ! T2 polygon: C00 → [e3 interior] → C10 → [e2 interior] → C11  (fan from C00)

        C00 = self % ctrlPts(P, 1, 1, :)
        C01 = self % ctrlPts(P, 1, 4, :)
        C10 = self % ctrlPts(P, 4, 1, :)
        C11 = self % ctrlPts(P, 4, 4, :)

        ! Collect sorted vertices for each edge (endpoints always included)
        P_adj = self % adj(P, 1)
        if (P_adj > 0 .and. wasSubdivided(P_adj)) then
          eP_idx = self % adjEdge(P, 1)
          call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                C00, C01, e1buf, n1)
        else
          n1 = 2;  e1buf(1,:) = C00;  e1buf(2,:) = C01
        end if

        P_adj = self % adj(P, 2)
        if (P_adj > 0 .and. wasSubdivided(P_adj)) then
          eP_idx = self % adjEdge(P, 2)
          call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                C10, C11, e2buf, n2)
        else
          n2 = 2;  e2buf(1,:) = C10;  e2buf(2,:) = C11
        end if

        P_adj = self % adj(P, 3)
        if (P_adj > 0 .and. wasSubdivided(P_adj)) then
          eP_idx = self % adjEdge(P, 3)
          call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                C00, C10, e3buf, n3)
        else
          n3 = 2;  e3buf(1,:) = C00;  e3buf(2,:) = C10
        end if

        P_adj = self % adj(P, 4)
        if (P_adj > 0 .and. wasSubdivided(P_adj)) then
          eP_idx = self % adjEdge(P, 4)
          call collectEdgeVerts(curPts, nCur, origPatch, uvRange, P_adj, eP_idx, &
                                C01, C11, e4buf, n4)
        else
          n4 = 2;  e4buf(1,:) = C01;  e4buf(2,:) = C11
        end if

        ! --- T1: polygon [e1_sorted || e4_sorted[2:]], fan from C00 ---
        nPoly = n1
        poly(1:n1, :) = e1buf(1:n1, :)
        do si = 2, n4
          nPoly = nPoly + 1
          poly(nPoly, :) = e4buf(si, :)
        end do
        do si = 2, nPoly - 1
          nTri = nTri + 1
          tris(nTri, 1, :) = poly(1, :)
          tris(nTri, 2, :) = poly(si, :)
          tris(nTri, 3, :) = poly(si+1, :)
        end do

        ! --- T2: polygon [e3_sorted || e2_sorted[2:]], fan from C00 ---
        nPoly = n3
        poly(1:n3, :) = e3buf(1:n3, :)
        do si = 2, n2
          nPoly = nPoly + 1
          poly(nPoly, :) = e2buf(si, :)
        end do
        do si = 2, nPoly - 1
          nTri = nTri + 1
          tris(nTri, 1, :) = poly(1, :)
          tris(nTri, 2, :) = poly(si, :)
          tris(nTri, 3, :) = poly(si+1, :)
        end do

      end if
    end do

    ! --- Step 5: Ray cast against the watertight triangle mesh ---
    hs = .not. self % rayCast(tris, nTri, r)

    deallocate(curPts, curWeights, origPatch, uvRange, tris)

    999 continue

  end function halfspace

  ! ---------------------------------------------------------------------------
  ! GJK convex hull containment test
  ! ---------------------------------------------------------------------------

  pure function cross3(a, b) result(c)
    real(defReal), dimension(3), intent(in) :: a, b
    real(defReal), dimension(3)             :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  pure function tripleCross(a, b, c) result(r)
    real(defReal), dimension(3), intent(in) :: a, b, c
    real(defReal), dimension(3)             :: r
    r = cross3(cross3(a, b), c)
  end function tripleCross

  !!
  !! GJK test: is point r inside the convex hull of pts(1:nPts,:)?
  !!
  pure function pointInConvexHull(pts, nPts, r) result(inside)
    real(defReal), dimension(:,:), intent(in) :: pts
    integer(shortInt), intent(in)             :: nPts
    real(defReal), dimension(3), intent(in)   :: r
    logical(defBool)                          :: inside
    real(defReal), dimension(4, 3) :: S
    integer(shortInt)              :: nS, i, bestIdx, iter
    real(defReal), dimension(3)    :: d, sup
    real(defReal)                  :: maxDot, dp
    integer(shortInt), parameter   :: MAX_ITER = 64

    inside = .false.
    if (nPts < 1) return

    ! Initial direction: from r toward centroid of the point set
    d = ZERO
    do i = 1, nPts
      d = d + pts(i,:)
    end do
    d = d / real(nPts, defReal) - r

    if (dot_product(d, d) < GJK_EPS) then
      inside = .true.
      return
    end if

    ! First support point
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

    if (dot_product(sup, d) < ZERO) return

    nS = 1
    S(1,:) = sup
    d = -sup

    do iter = 1, MAX_ITER
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

      if (dot_product(sup, d) < ZERO) return

      nS = nS + 1
      S(nS,:) = sup

      call gjkDoSimplex(S, nS, d, inside)
      if (inside) return
    end do

    inside = .true.  ! failed to converge: conservatively inside

  end function pointInConvexHull

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

    if (nS == 4) then
      a  = S(4,:);  b = S(3,:);  c = S(2,:);  dd = S(1,:)
      ab = b - a;   ac = c - a;  ad = dd - a;  ao = -a

      abc = cross3(ab, ac);  if (dot_product(abc, ad) > ZERO) abc = -abc
      acd = cross3(ac, ad);  if (dot_product(acd, ab) > ZERO) acd = -acd
      adb = cross3(ad, ab);  if (dot_product(adb, ac) > ZERO) adb = -adb

      if (dot_product(abc, ao) > ZERO) then
        if (dot_product(cross3(ab, ac), ad) < ZERO) then
          S(1,:) = c;  S(2,:) = b;  S(3,:) = a
        else
          S(1,:) = b;  S(2,:) = c;  S(3,:) = a
        end if
        nS = 3
      else if (dot_product(acd, ao) > ZERO) then
        if (dot_product(cross3(ac, ad), ab) < ZERO) then
          S(1,:) = dd;  S(2,:) = c;  S(3,:) = a
        else
          S(1,:) = c;  S(2,:) = dd;  S(3,:) = a
        end if
        nS = 3
      else if (dot_product(adb, ao) > ZERO) then
        if (dot_product(cross3(ad, ab), ac) < ZERO) then
          S(1,:) = b;  S(2,:) = dd;  S(3,:) = a
        else
          S(1,:) = dd;  S(2,:) = b;  S(3,:) = a
        end if
        nS = 3
      else
        inside = .true.
        return
      end if
    end if

    if (nS == 3) then
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
      return
    end if

    if (nS == 2) then
      a  = S(2,:);  b = S(1,:)
      ab = b - a;   ao = -a

      if (dot_product(ab, ao) > ZERO) then
        tcross = tripleCross(ab, ao, ab)
        if (dot_product(tcross, tcross) < GJK_EPS) then
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
      return
    end if

  end subroutine gjkDoSimplex

  ! ---------------------------------------------------------------------------
  ! Bezier patch subdivision
  ! ---------------------------------------------------------------------------

  !!
  !! Rational Bezier patch subdivision at t=0.5 using homogeneous De Casteljau.
  !! Lifts control points to 4D homogeneous space (wx,wy,wz,w), subdivides with
  !! standard linear averaging (correct for any dimensionality), then dehomogenizes.
  !! For unit weights this is identical to polynomial De Casteljau.
  !!
  subroutine subdividePatch(patch, wts, Q00, Q01, Q10, Q11, Q00w, Q01w, Q10w, Q11w)
    real(defReal), dimension(4,4,3), intent(in)  :: patch
    real(defReal), dimension(4,4),   intent(in)  :: wts
    real(defReal), dimension(4,4,3), intent(out) :: Q00, Q01, Q10, Q11
    real(defReal), dimension(4,4),   intent(out) :: Q00w, Q01w, Q10w, Q11w
    real(defReal), dimension(4,4,4) :: H, leftH, rightH, HL, HR
    real(defReal), dimension(4,4)   :: rL, rR
    integer(shortInt) :: i, j

    ! Lift to 4D homogeneous: H(i,j,1:3) = w*P, H(i,j,4) = w
    do i = 1, 4
      do j = 1, 4
        H(i,j,1:3) = wts(i,j) * patch(i,j,:)
        H(i,j,4)   = wts(i,j)
      end do
    end do

    ! Subdivide in v-direction (each u-row H(i,:,:))
    do i = 1, 4
      call subdivideRowH(H(i,:,:), rL, rR)
      leftH(i,:,:)  = rL
      rightH(i,:,:) = rR
    end do

    ! Subdivide left-v half in u-direction (each v-column leftH(:,j,:))
    do j = 1, 4
      call subdivideRowH(leftH(:,j,:), rL, rR)
      HL(:,j,:) = rL   ! Q00 homogeneous
      HR(:,j,:) = rR   ! Q10 homogeneous
    end do
    call dehomogenize4x4(HL, Q00, Q00w)
    call dehomogenize4x4(HR, Q10, Q10w)

    ! Subdivide right-v half in u-direction
    do j = 1, 4
      call subdivideRowH(rightH(:,j,:), rL, rR)
      HL(:,j,:) = rL   ! Q01 homogeneous
      HR(:,j,:) = rR   ! Q11 homogeneous
    end do
    call dehomogenize4x4(HL, Q01, Q01w)
    call dehomogenize4x4(HR, Q11, Q11w)

  end subroutine subdividePatch

  !! De Casteljau subdivision at t=0.5 for a row of 4 homogeneous points (4-component).
  subroutine subdivideRowH(row, left, right)
    real(defReal), dimension(4,4), intent(in)  :: row
    real(defReal), dimension(4,4), intent(out) :: left, right
    real(defReal), dimension(4) :: p01, p12, p23, p012, p123, p0123

    p01   = HALF * (row(1,:) + row(2,:))
    p12   = HALF * (row(2,:) + row(3,:))
    p23   = HALF * (row(3,:) + row(4,:))
    p012  = HALF * (p01  + p12)
    p123  = HALF * (p12  + p23)
    p0123 = HALF * (p012 + p123)

    left(1,:)  = row(1,:);  left(2,:)  = p01;    left(3,:)  = p012;  left(4,:)  = p0123
    right(1,:) = p0123;     right(2,:) = p123;   right(3,:) = p23;   right(4,:) = row(4,:)

  end subroutine subdivideRowH

  !! Dehomogenize a 4x4 array of 4D homogeneous points to 3D positions and weights.
  subroutine dehomogenize4x4(H, pts, wts)
    real(defReal), dimension(4,4,4), intent(in)  :: H
    real(defReal), dimension(4,4,3), intent(out) :: pts
    real(defReal), dimension(4,4),   intent(out) :: wts
    integer(shortInt) :: i, j
    real(defReal), parameter :: W_MIN = 1.0E-30_defReal

    do i = 1, 4
      do j = 1, 4
        wts(i,j) = H(i,j,4)
        if (abs(wts(i,j)) > W_MIN) then
          pts(i,j,:) = H(i,j,1:3) / wts(i,j)
        else
          pts(i,j,:) = ZERO
        end if
      end do
    end do

  end subroutine dehomogenize4x4

  ! ---------------------------------------------------------------------------
  ! Ray casting
  ! ---------------------------------------------------------------------------

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
    real(defReal), parameter    :: RAY_NUDGE = 1.0E-5_defReal

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
      rLocal = r + (/ RAY_NUDGE, RAY_NUDGE, RAY_NUDGE /)
      count = 0
      do i = 1, nTri
        v0 = tris(i, 1, :);  v1 = tris(i, 2, :);  v2 = tris(i, 3, :)
        call self % rayTriangle(rLocal, rayDir, v0, v1, v2, hit, nearZero)
        if (hit) count = count + 1
      end do
    end if

    inside = mod(count, 2) == 1

  end function rayCast

  subroutine rayCastDebug(self, tris, nTri, r)
    class(bezierVolume), intent(in)              :: self
    real(defReal), dimension(:,:,:), intent(in)  :: tris
    integer(shortInt), intent(in)                :: nTri
    real(defReal), dimension(3), intent(in)      :: r
    real(defReal), dimension(3) :: rayDir, v0, v1, v2
    integer(shortInt)           :: count, i
    logical(defBool)            :: hit, nz

    rayDir = (/ ONE, ONE / 3.0_defReal, ONE / 7.0_defReal /)
    rayDir = rayDir / norm2(rayDir)

    count = 0
    write(*, '(A, 3F10.6)') '  Ray direction: ', rayDir
    write(*, '(A, I8)')     '  Total triangles: ', nTri

    do i = 1, nTri
      v0 = tris(i, 1, :);  v1 = tris(i, 2, :);  v2 = tris(i, 3, :)
      call self % rayTriangle(r, rayDir, v0, v1, v2, hit, nz)
      if (hit) then
        count = count + 1
        write(*, '(A, I6)') '  HIT tri: ', i
        write(*, '(A, 3F8.4)') '    v0: ', v0
        write(*, '(A, 3F8.4)') '    v1: ', v1
        write(*, '(A, 3F8.4)') '    v2: ', v2
      end if
    end do

    write(*, '(A, I4, A, L3)') '  Total hits: ', count, '  inside=', mod(count, 2) == 1

  end subroutine rayCastDebug

  !!
  !! Watertight ray-triangle intersection (Woop, Benthin & Wald 2013).
  !! Edge functions depend only on each edge's two vertices, so shared edges
  !! between adjacent triangles produce bitwise-identical results (no gaps).
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
    real(defReal), parameter :: EPS           = 1.0E-10_defReal
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

    if (e0 < -EPS .or. e1 < -EPS .or. e2 < -EPS) then
      if (e0 > EPS .or. e1 > EPS .or. e2 > EPS) return
    end if

    det = e0 + e1 + e2
    if (abs(det) < EPS) then
      nearZero = .true.
      return
    end if

    Az = Sz * A(kz);  Bz = Sz * B(kz);  Cz = Sz * C(kz)
    t  = e0 * Az + e1 * Bz + e2 * Cz

    if (abs(t) < NEAR_ZERO_TOL * abs(det)) nearZero = .true.

    if (det > ZERO) then
      if (t < EPS * det) return
    else
      if (t > EPS * det) return
    end if

    hit = .true.

  end subroutine rayTriangle

  ! ---------------------------------------------------------------------------
  ! Module-level helper functions
  ! ---------------------------------------------------------------------------

  !!
  !! Collect all unique corner points from sub-patches of adjPatch that lie on
  !! edge edgeOnAdj of the original patch, sorted by distance from ptA toward ptB.
  !! Always includes ptA (first) and ptB (last) — buf has at least 2 entries.
  !!
  !!
  !! Collect T-junction vertices on one edge of sub-patch Si from finer sub-patches
  !! of the same original patch (within-patch T-junction detection).
  !!
  !! eP_Si edge convention matches the patch edge numbering:
  !!   1 = u=uMin (left),  2 = u=uMax (right)
  !!   3 = v=vMin (bottom), 4 = v=vMax (top)
  !!
  !! ptA and ptB are the 3D end-points of Si's edge in the direction the polygon
  !! fan expects them (matches the T1/T2 polygon winding in halfspace).
  !!
  !! Returns buf sorted ptA→ptB with all T-junction vertices inserted.
  !! Always contains at least ptA and ptB.
  !!
  subroutine collectIntraPatchEdgeVerts(curPts, nCur, origPatch, uvRange, &
                                         Si_idx, eP_Si, ptA, ptB, buf, nBuf)
    real(defReal), dimension(:,:,:,:), intent(in) :: curPts
    integer(shortInt), intent(in)                 :: nCur
    integer(shortInt), dimension(:), intent(in)   :: origPatch
    real(defReal), dimension(:,:), intent(in)     :: uvRange
    integer(shortInt), intent(in)                 :: Si_idx, eP_Si
    real(defReal), dimension(3), intent(in)       :: ptA, ptB
    real(defReal), dimension(:,:), intent(out)    :: buf
    integer(shortInt), intent(out)                :: nBuf

    integer(shortInt) :: i, si, P
    real(defReal), dimension(3)       :: pt1, pt2, edgeDir, tmpV
    real(defReal), dimension(size(buf,1)) :: dotArr
    real(defReal) :: tmpD, boundaryVal, pMin, pMax
    logical(defBool) :: onEdge, isDup
    real(defReal), parameter :: EV_TOL = 1.0E-10_defReal

    P       = origPatch(Si_idx)
    nBuf    = 0
    edgeDir = ptB - ptA

    ! UV boundary value and range for this edge of Si
    select case (eP_Si)
      case(1);  boundaryVal = uvRange(Si_idx,1);  pMin = uvRange(Si_idx,3);  pMax = uvRange(Si_idx,4)
      case(2);  boundaryVal = uvRange(Si_idx,2);  pMin = uvRange(Si_idx,3);  pMax = uvRange(Si_idx,4)
      case(3);  boundaryVal = uvRange(Si_idx,3);  pMin = uvRange(Si_idx,1);  pMax = uvRange(Si_idx,2)
      case(4);  boundaryVal = uvRange(Si_idx,4);  pMin = uvRange(Si_idx,1);  pMax = uvRange(Si_idx,2)
    end select

    do i = 1, nCur
      if (i == Si_idx) cycle
      if (origPatch(i) /= P) cycle

      ! Check that Sj shares Si's edge boundary and lies within Si's range
      onEdge = .false.
      select case (eP_Si)
        case(1)
          onEdge = (abs(uvRange(i,2) - boundaryVal) < EV_TOL) .and. &
                   (uvRange(i,3) >= pMin - EV_TOL) .and. &
                   (uvRange(i,4) <= pMax + EV_TOL)
        case(2)
          onEdge = (abs(uvRange(i,1) - boundaryVal) < EV_TOL) .and. &
                   (uvRange(i,3) >= pMin - EV_TOL) .and. &
                   (uvRange(i,4) <= pMax + EV_TOL)
        case(3)
          onEdge = (abs(uvRange(i,4) - boundaryVal) < EV_TOL) .and. &
                   (uvRange(i,1) >= pMin - EV_TOL) .and. &
                   (uvRange(i,2) <= pMax + EV_TOL)
        case(4)
          onEdge = (abs(uvRange(i,3) - boundaryVal) < EV_TOL) .and. &
                   (uvRange(i,1) >= pMin - EV_TOL) .and. &
                   (uvRange(i,2) <= pMax + EV_TOL)
      end select
      if (.not. onEdge) cycle

      ! Extract the 2 corners of Sj that lie on the shared boundary
      ! Si's e=1 (left) is adjacent to Sj's e=2 (right), etc.
      select case (eP_Si)
        case(1);  pt1 = curPts(i,4,1,:);  pt2 = curPts(i,4,4,:)
        case(2);  pt1 = curPts(i,1,1,:);  pt2 = curPts(i,1,4,:)
        case(3);  pt1 = curPts(i,1,4,:);  pt2 = curPts(i,4,4,:)
        case(4);  pt1 = curPts(i,1,1,:);  pt2 = curPts(i,4,1,:)
      end select

      isDup = .false.
      do si = 1, nBuf
        if (norm2(buf(si,:) - pt1) < EV_TOL) then;  isDup = .true.;  exit;  end if
      end do
      if (.not. isDup .and. nBuf < size(buf,1)) then
        nBuf = nBuf + 1;  buf(nBuf,:) = pt1
      end if

      isDup = .false.
      do si = 1, nBuf
        if (norm2(buf(si,:) - pt2) < EV_TOL) then;  isDup = .true.;  exit;  end if
      end do
      if (.not. isDup .and. nBuf < size(buf,1)) then
        nBuf = nBuf + 1;  buf(nBuf,:) = pt2
      end if
    end do

    if (nBuf < 2) then
      nBuf = 2;  buf(1,:) = ptA;  buf(2,:) = ptB
      return
    end if

    ! Insertion sort by dot product along edge direction ptA→ptB
    do si = 1, nBuf
      dotArr(si) = dot_product(buf(si,:) - ptA, edgeDir)
    end do
    do si = 2, nBuf
      tmpV = buf(si,:);  tmpD = dotArr(si)
      i = si - 1
      do while (i >= 1 .and. dotArr(i) > tmpD)
        buf(i+1,:) = buf(i,:);  dotArr(i+1) = dotArr(i)
        i = i - 1
      end do
      buf(i+1,:) = tmpV;  dotArr(i+1) = tmpD
    end do

  end subroutine collectIntraPatchEdgeVerts

  subroutine collectEdgeVerts(curPts, nCur, origPatch, uvRange, adjPatch, edgeOnAdj, &
                              ptA, ptB, buf, nBuf)
    real(defReal), dimension(:,:,:,:), intent(in) :: curPts
    integer(shortInt), intent(in)                 :: nCur
    integer(shortInt), dimension(:), intent(in)   :: origPatch
    real(defReal), dimension(:,:), intent(in)     :: uvRange
    integer(shortInt), intent(in)                 :: adjPatch, edgeOnAdj
    real(defReal), dimension(3), intent(in)       :: ptA, ptB
    real(defReal), dimension(:,:), intent(out)    :: buf
    integer(shortInt), intent(out)                :: nBuf

    integer(shortInt) :: i, si
    logical(defBool)  :: onEdge, isDup
    real(defReal), dimension(3) :: pt1, pt2, edgeDir, tmpV
    real(defReal), dimension(size(buf,1)) :: dotArr
    real(defReal) :: tmpD
    real(defReal), parameter :: EV_TOL = 1.0E-10_defReal

    nBuf = 0
    edgeDir = ptB - ptA

    do i = 1, nCur
      if (origPatch(i) /= adjPatch) cycle

      onEdge = .false.
      select case (edgeOnAdj)
        case(1);  onEdge = (uvRange(i,1) < EV_TOL)
        case(2);  onEdge = (uvRange(i,2) > ONE - EV_TOL)
        case(3);  onEdge = (uvRange(i,3) < EV_TOL)
        case(4);  onEdge = (uvRange(i,4) > ONE - EV_TOL)
      end select
      if (.not. onEdge) cycle

      select case (edgeOnAdj)
        case(1);  pt1 = curPts(i,1,1,:);  pt2 = curPts(i,1,4,:)
        case(2);  pt1 = curPts(i,4,1,:);  pt2 = curPts(i,4,4,:)
        case(3);  pt1 = curPts(i,1,1,:);  pt2 = curPts(i,4,1,:)
        case(4);  pt1 = curPts(i,1,4,:);  pt2 = curPts(i,4,4,:)
      end select

      isDup = .false.
      do si = 1, nBuf
        if (norm2(buf(si,:) - pt1) < EV_TOL) then;  isDup = .true.;  exit;  end if
      end do
      if (.not. isDup .and. nBuf < size(buf,1)) then
        nBuf = nBuf + 1;  buf(nBuf,:) = pt1
      end if

      isDup = .false.
      do si = 1, nBuf
        if (norm2(buf(si,:) - pt2) < EV_TOL) then;  isDup = .true.;  exit;  end if
      end do
      if (.not. isDup .and. nBuf < size(buf,1)) then
        nBuf = nBuf + 1;  buf(nBuf,:) = pt2
      end if
    end do

    ! Fallback: if no sub-patches found, just return the two endpoints
    if (nBuf < 2) then
      nBuf = 2;  buf(1,:) = ptA;  buf(2,:) = ptB
      return
    end if

    ! Insertion sort by dot product with edge direction (ptA → ptB)
    do si = 1, nBuf
      dotArr(si) = dot_product(buf(si,:) - ptA, edgeDir)
    end do
    do si = 2, nBuf
      tmpV = buf(si,:);  tmpD = dotArr(si)
      i = si - 1
      do while (i >= 1 .and. dotArr(i) > tmpD)
        buf(i+1,:) = buf(i,:);  dotArr(i+1) = dotArr(i)
        i = i - 1
      end do
      buf(i+1,:) = tmpV;  dotArr(i+1) = tmpD
    end do

  end subroutine collectEdgeVerts

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
    real(defReal), parameter :: EV_TOL = 1.0E-10_defReal

    onBoundary = .false.
    select case (edge)
      case(1);  onBoundary = uvRange(idx,1) < EV_TOL
      case(2);  onBoundary = uvRange(idx,2) > ONE - EV_TOL
      case(3);  onBoundary = uvRange(idx,3) < EV_TOL
      case(4);  onBoundary = uvRange(idx,4) > ONE - EV_TOL
    end select

  end function edgeOnOuterBoundary

  !!
  !! Checks whether an original patch is eligible for the AABB-only subdivision fallback
  !!
  !! Scopes the Step 3 AABB-only fallback (triggered when the strict hull test
  !! finds no patch, typically floating-point precision on an already-tight,
  !! deeply-subdivided hull) so it does not pull in unrelated patches purely by
  !! loose bounding-box coincidence. See the Step 3 call site for the full
  !! rationale and the Gumbo case that motivated it.
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
  !! Merges a neighbouring sub-patch's boundary vertices into this sub-patch's buffer
  !!
  !! Reconciles two adjacent patches' independently-generated T-junction vertices
  !! along a shared edge when both were selectively subdivided (see call site for
  !! why this is needed -- collectIntraPatchEdgeVerts alone only sees one side).
  !! Vertices from extraBuf are accepted only if their projection onto ptA->ptB
  !! falls within [0,1] (extraBuf may cover more of the shared edge than this
  !! sub-patch's own [ptA,ptB] portion of it); accepted vertices are merged into
  !! buf with exact-position deduplication, then the whole of buf is re-sorted
  !! by position along ptA->ptB.
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
    logical(defBool)  :: isDup
    real(defReal), dimension(3) :: edgeDir, tmpV
    real(defReal), dimension(size(buf,1)) :: dotArr
    real(defReal) :: tmpD, tParam, lenSq
    real(defReal), parameter :: EV_TOL = 1.0E-10_defReal
    real(defReal), parameter :: RANGE_TOL = 1.0E-6_defReal

    edgeDir = ptB - ptA
    lenSq = dot_product(edgeDir, edgeDir)

    do sj = 1, nExtra
      ! extraBuf may come from a neighbour patch's FULL boundary edge, which can
      ! span more than THIS sub-patch's own [ptA,ptB] portion of it (e.g. when
      ! the neighbour is subdivided into multiple pieces along the shared edge,
      ! each call here only owns one piece). Reject anything whose projection
      ! falls outside the ptA-ptB span. No perpendicular-distance check here:
      ! this is a genuinely CURVED shared boundary, so interior points do not
      ! lie on the straight chord between ptA and ptB -- every candidate from
      ! collectEdgeVerts is already known to sit exactly on the true shared
      ! curve (buildAdjacency verified the two patches' boundary control points
      ! match exactly), so projection range is the only thing worth checking.
      if (lenSq > EV_TOL) then
        tParam = dot_product(extraBuf(sj,:) - ptA, edgeDir) / lenSq
        if (tParam < -RANGE_TOL .or. tParam > ONE + RANGE_TOL) cycle
      end if

      isDup = .false.
      do si = 1, nBuf
        if (norm2(buf(si,:) - extraBuf(sj,:)) < EV_TOL) then
          isDup = .true.
          exit
        end if
      end do
      if (.not. isDup .and. nBuf < size(buf,1)) then
        nBuf = nBuf + 1
        buf(nBuf,:) = extraBuf(sj,:)
      end if
    end do

    ! Re-sort by position along ptA->ptB
    do si = 1, nBuf
      dotArr(si) = dot_product(buf(si,:) - ptA, edgeDir)
    end do
    do si = 2, nBuf
      tmpV = buf(si,:);  tmpD = dotArr(si)
      sj = si - 1
      do while (sj >= 1 .and. dotArr(sj) > tmpD)
        buf(sj+1,:) = buf(sj,:);  dotArr(sj+1) = dotArr(sj)
        sj = sj - 1
      end do
      buf(sj+1,:) = tmpV;  dotArr(sj+1) = tmpD
    end do

  end subroutine mergeEdgeVerts

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

  !!
  !! Compute the De Casteljau midpoint (t=0.5) of a patch's boundary edge.
  !!
  pure function patchEdgeMid(patch, e) result(M)
    real(defReal), dimension(4,4,3), intent(in) :: patch
    integer(shortInt), intent(in)               :: e
    real(defReal), dimension(3)                 :: M
    real(defReal), dimension(4,3)               :: edge
    real(defReal), dimension(3)                 :: p01, p12, p23, p012, p123

    edge = patchEdge(patch, e)

    p01  = HALF * (edge(1,:) + edge(2,:))
    p12  = HALF * (edge(2,:) + edge(3,:))
    p23  = HALF * (edge(3,:) + edge(4,:))
    p012 = HALF * (p01 + p12)
    p123 = HALF * (p12 + p23)
    M    = HALF * (p012 + p123)

  end function patchEdgeMid

  ! ---------------------------------------------------------------------------
  ! Kill
  ! ---------------------------------------------------------------------------

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

end module bezierVolume_class
