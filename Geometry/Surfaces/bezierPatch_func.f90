module bezierPatch_func
  !! Pure evaluation/subdivision routines for a single (possibly rational)
  !! bicubic Bezier patch: a 4x4 grid of 3D control points, each with an
  !! optional rational weight (weight 1 everywhere reduces to a plain
  !! polynomial Bezier patch).
  !!
  !! All routines work via De Casteljau's algorithm, lifted into 4D
  !! homogeneous space (homogeneousLift/dehomogenize) so that rational
  !! (NURBS-style) weights are handled by the same lerp-based recursion as
  !! the unweighted case, with no separate rational code path.
  !!
  !! Deliberately a free-standing, stateless module rather than private
  !! procedures inside a class: two independent callers need this same
  !! patch-level math but have no shared class instance to hang it off --
  !! bezierVolume_class (Geometry/Surfaces/bezierVolume_class.f90) calls
  !! subdividePatch while adaptively subdividing patches for its halfspace
  !! test, and Apps/sdfSampler.f90 calls evalPatch directly on raw
  !! control points/weights parsed from its own input dict (it never
  !! constructs a bezierVolume, so it cannot reach this code as a type-bound
  !! procedure) to place near-surface training points at exact (u,v)
  !! parametric locations for its banded sampling strategy.
  !!
  !! Public Members:
  !!   homogeneousLift  -> lift a patch's control points + weights to 4D
  !!   dehomogenize     -> inverse of homogeneousLift
  !!   deCasteljauLerp  -> evaluate one row of 4 homogeneous points at t
  !!   deCasteljauSplit -> split one row of 4 homogeneous points at t into
  !!                       its two sub-curves' control points
  !!   evalPatch        -> evaluate the full patch at parameter (u,v)
  !!   subdividePatch   -> split the full patch into its 4 sub-patches at
  !!                       (u,v)
  !!

  use numPrecision
  implicit none
  private

  public :: homogeneousLift
  public :: dehomogenize
  public :: deCasteljauLerp
  public :: deCasteljauSplit
  public :: evalPatch
  public :: subdividePatch

  real(defReal), parameter :: MIN_WEIGHT = 1.0E-30_defReal

contains

  !!
  !! Lift a bicubic patch's control points and rational weights to 4D
  !! homogeneous space: H(i,j,1:3) = w*P, H(i,j,4) = w.
  !!
  pure function homogeneousLift(ctrlPts, wts) result(H)
    real(defReal), dimension(4,4,3), intent(in) :: ctrlPts
    real(defReal), dimension(4,4), intent(in)   :: wts
    real(defReal), dimension(4,4,4) :: H
    integer(shortInt) :: i, j

    do i = 1, 4
      do j = 1, 4
        H(i,j,1:3) = wts(i,j) * ctrlPts(i,j,:)
        H(i,j,4)   = wts(i,j)
      end do
    end do

  end function homogeneousLift

  !!
  !! Inverse of homogeneousLift: recover 3D positions and rational weights
  !! from a 4x4 array of 4D homogeneous points.
  !!
  subroutine dehomogenize(H, ctrlPts, wts)
    real(defReal), dimension(4,4,4), intent(in)  :: H
    real(defReal), dimension(4,4,3), intent(out) :: ctrlPts
    real(defReal), dimension(4,4),   intent(out) :: wts
    integer(shortInt) :: i, j

    do i = 1, 4
      do j = 1, 4
        wts(i,j) = H(i,j,4)
        if (abs(wts(i,j)) > MIN_WEIGHT) then
          ctrlPts(i,j,:) = H(i,j,1:3) / wts(i,j)
        else
          ctrlPts(i,j,:) = ZERO
        end if
      end do
    end do

  end subroutine dehomogenize

  !!
  !! Linear-interpolation reduction of 4 homogeneous points to a single point
  !! at parameter t, via 3 iterative lerps (cubic De Casteljau).
  !!
  !! Args:
  !!   row [in] -> 4 control points, each with 4 homogeneous components
  !!   t [in] -> Parameter in [0,1] at which to evaluate
  !!
  !! Result:
  !!   The 4-component homogeneous point at parameter t
  !!
  pure function deCasteljauLerp(row, t) result(q)
    real(defReal), dimension(4,4), intent(in) :: row
    real(defReal), intent(in)                 :: t
    real(defReal), dimension(4)               :: q
    real(defReal), dimension(4) :: p01, p12, p23, p012, p123

    p01  = (ONE - t) * row(1,:) + t * row(2,:)
    p12  = (ONE - t) * row(2,:) + t * row(3,:)
    p23  = (ONE - t) * row(3,:) + t * row(4,:)
    p012 = (ONE - t) * p01  + t * p12
    p123 = (ONE - t) * p12  + t * p23
    q    = (ONE - t) * p012 + t * p123

  end function deCasteljauLerp

  !!
  !! De Casteljau subdivision of a row of 4 homogeneous points at parameter t,
  !! into the left and right control-point sequences of the two resulting
  !! sub-curves. At t=0.5 this reduces to standard midpoint subdivision.
  !!
  !! Args:
  !!   row [in] -> 4 control points, each with 4 homogeneous components
  !!   t [in] -> Split parameter in [0,1]
  !!
  !! Result:
  !!   left [out], right [out] -> the two resulting 4-point control sequences
  !!
  subroutine deCasteljauSplit(row, t, left, right)
    real(defReal), dimension(4,4), intent(in)  :: row
    real(defReal), intent(in)                  :: t
    real(defReal), dimension(4,4), intent(out) :: left, right
    real(defReal), dimension(4) :: p01, p12, p23, p012, p123, p0123

    p01   = (ONE - t) * row(1,:) + t * row(2,:)
    p12   = (ONE - t) * row(2,:) + t * row(3,:)
    p23   = (ONE - t) * row(3,:) + t * row(4,:)
    p012  = (ONE - t) * p01  + t * p12
    p123  = (ONE - t) * p12  + t * p23
    p0123 = (ONE - t) * p012 + t * p123

    left(1,:)  = row(1,:);  left(2,:)  = p01;    left(3,:)  = p012;  left(4,:)  = p0123
    right(1,:) = p0123;     right(2,:) = p123;   right(3,:) = p23;   right(4,:) = row(4,:)

  end subroutine deCasteljauSplit

  !!
  !! Evaluate a (rational) bicubic Bezier patch at parameter (uParam, vParam)
  !! via homogeneous-lift De Casteljau.
  !!
  !! Args:
  !!   ctrlPts [in] -> Patch control points, ctrlPts(i,j,:); i = u-index, j = v-index
  !!   wts [in] -> Patch rational weights, wts(i,j)
  !!   uParam [in] -> Parameter in [0,1] along the u-direction
  !!   vParam [in] -> Parameter in [0,1] along the v-direction
  !!
  !! Result:
  !!   The 3D surface point at (uParam, vParam)
  !!
  pure function evalPatch(ctrlPts, wts, uParam, vParam) result(p)
    real(defReal), dimension(4,4,3), intent(in) :: ctrlPts
    real(defReal), dimension(4,4), intent(in)   :: wts
    real(defReal), intent(in)                   :: uParam, vParam
    real(defReal), dimension(3)                 :: p
    real(defReal), dimension(4,4,4) :: H
    real(defReal), dimension(4,4)   :: Hu
    real(defReal), dimension(4)     :: Hp
    integer(shortInt) :: i

    H = homogeneousLift(ctrlPts, wts)

    do i = 1, 4
      Hu(i,:) = deCasteljauLerp(H(i,:,:), vParam)
    end do

    Hp = deCasteljauLerp(Hu, uParam)

    if (abs(Hp(4)) > MIN_WEIGHT) then
      p = Hp(1:3) / Hp(4)
    else
      p = ZERO
    end if

  end function evalPatch

  !!
  !! Subdivide a (rational) bicubic Bezier patch into its 4 sub-patches at
  !! split parameter (uParam, vParam), via homogeneous-lift De Casteljau.
  !!
  !! Args:
  !!   ctrlPts [in], wts [in] -> Patch control points and rational weights
  !!   uParam [in], vParam [in] -> Split parameters in [0,1]
  !!
  !! Result:
  !!   Q00,Q01,Q10,Q11 [out] / Q00w,Q01w,Q10w,Q11w [out] -> the 4 sub-patches'
  !!     control points and weights, indexed [uMin/uMax][vMin/vMax]
  !!
  subroutine subdividePatch(ctrlPts, wts, uParam, vParam, &
                            Q00, Q01, Q10, Q11, Q00w, Q01w, Q10w, Q11w)
    real(defReal), dimension(4,4,3), intent(in)  :: ctrlPts
    real(defReal), dimension(4,4),   intent(in)  :: wts
    real(defReal), intent(in)                    :: uParam, vParam
    real(defReal), dimension(4,4,3), intent(out) :: Q00, Q01, Q10, Q11
    real(defReal), dimension(4,4),   intent(out) :: Q00w, Q01w, Q10w, Q11w
    real(defReal), dimension(4,4,4) :: H, leftH, rightH, HL, HR
    real(defReal), dimension(4,4)   :: rL, rR
    integer(shortInt) :: i, j

    H = homogeneousLift(ctrlPts, wts)

    ! Subdivide in v-direction (each u-row H(i,:,:))
    do i = 1, 4
      call deCasteljauSplit(H(i,:,:), vParam, rL, rR)
      leftH(i,:,:)  = rL
      rightH(i,:,:) = rR
    end do

    ! Subdivide left-v half in u-direction (each v-column leftH(:,j,:))
    do j = 1, 4
      call deCasteljauSplit(leftH(:,j,:), uParam, rL, rR)
      HL(:,j,:) = rL   ! Q00 homogeneous
      HR(:,j,:) = rR   ! Q10 homogeneous
    end do
    call dehomogenize(HL, Q00, Q00w)
    call dehomogenize(HR, Q10, Q10w)

    ! Subdivide right-v half in u-direction
    do j = 1, 4
      call deCasteljauSplit(rightH(:,j,:), uParam, rL, rR)
      HL(:,j,:) = rL   ! Q01 homogeneous
      HR(:,j,:) = rR   ! Q11 homogeneous
    end do
    call dehomogenize(HL, Q01, Q01w)
    call dehomogenize(HR, Q11, Q11w)

  end subroutine subdividePatch

end module bezierPatch_func
