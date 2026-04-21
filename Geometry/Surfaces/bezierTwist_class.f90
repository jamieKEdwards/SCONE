module bezierTwist_class

  use numPrecision
  use universalVariables, only : X_AXIS, Y_AXIS, Z_AXIS, INF, SURF_TOL
  use genericProcedures,  only : fatalError, dotProduct, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill
  implicit none
  private

  !!
  !!
  !! Bezier shape defined by intersection watertight bezier curves with twisted extrusion in Z-direction
  !!
  !! In 2D, given a watertight set fo bezier curves, this surface has 
  !! halfspace, but not distance too. This is work in progress. 
  !! The class uses various point in polygon tests to see if a point lies 
  !! within  the control point boundaries of a curve or collection of curves. 
  !! This is then followed by a Graham Scan to ensure we have the convex 
  !! hull of each bezier curve, followed by a recursive subdivision to check
  !! which side of the curve a point lies. 
  !!
  !! For now the control points should be inputted in a list in order of ( curve((point)(point)) curve((point)(point)) )
  !! the order of curve should be provided such that scone can divide up the points correctly. 
  !!
  !! This function is identical to the bezier Shape function with the adition of a twisting transformation along the z axis
  !! The additional input to allow the twist is some twisting period (in cm) across which the shape will rotate by pi
  !!
  !! Surface tolerance: SURF_TOL
  !! Nudge: NUDGE
  !!
  !! Sample dictionary input:
  !!  pl { type bezierTwist; id 16; order 3; minZ -10; maxZ 10; twist 10; ctrlPts   (1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0,
  !!                                                    1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0,
  !!                                                    1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0);
  !!      }
  !!
  !! Private members:
  !!   norm -> Normal vector (normalised) [c1, c2, c3]
  !!   offset -> Offset, c4 (also normalised)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: bezierTwist
    private
    real(defReal), dimension(3)                     :: norm = ZERO
    real(defReal)                                   :: offset = ZERO
    real(defReal)                                   :: twist = ONE
    real(defReal)                                   :: minZ = -ONE
    real(defReal)                                   :: maxZ = ONE
    real(defReal), dimension(:, :, :), allocatable  :: ctrlPts
    integer(shortInt)                               :: numCurves = 0
    integer(shortInt)                               :: order = 0
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
    procedure :: subdivide
    procedure :: order3Hull
    procedure :: grahamScan
    procedure :: orientation
    procedure :: inOrOut
    procedure :: inConvexPoly
    procedure :: findInnerOuterPaths
  end type bezierTwist


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(bezierTwist), intent(in)  :: self
    character(:), allocatable  :: str

    str = 'bezierTwist'

  end function myType

  !!
  !! Initialise bezierTwist from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if id < 0.
  !!   fatalError if is not a bezierTwist (coeffcients 1-3 are 0.0)
  !!
  subroutine init(self, dict)
    class(bezierTwist), intent(inout)        :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id, n, m, i, j, k, order
    real(defReal), dimension(:), allocatable :: ctrlPtsList
    real(defReal)                            :: twist, minZ, maxZ
    character(100), parameter :: Here = 'init (bezierTwist_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(ctrlPtsList,'ctrlPts')
    call dict % get(order,'order')
    call dict % get(twist,'twist')
    call dict % get(maxZ,'maxZ')
    call dict % get(minZ,'minZ')
    
    n = size(ctrlPtsList)    
    self % order = order
    self % numCurves = (n / 3) / (self % order + 1)
    self % twist = twist
    self % minZ = minZ
    self % maxZ = maxZ

    ! Check values
    if (id < 1) then
      call fatalError(Here,'Invalid surface id provided. ID must be > 1')
    end if

    ! Load data
    call self % setID(id)
    
    ! Unpack points and fill array of correct dimensions
    allocate(self % ctrlPts(self % numCurves, self % order + 1, 3))

    m=1
    do while (m<=n) 
      do i=1, self % numCurves
        do j=1, self % order + 1
          do k=1, 3
            self % ctrlPts(i, j, k) = ctrlPtsList(m)
            m = m + 1
          end do 
        end do
      end do
    end do

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  !! Returns value of max x, y for the given points and inf in z 
  !! 
  pure function boundingBox(self) result(aabb)
    class(bezierTwist), intent(in)  :: self
    real(defReal), dimension(6)     :: aabb
    integer(shortInt)               :: i, j 

    aabb(1) = INF
    aabb(2) = INF
    aabb(3) = -INF
    aabb(4) = -INF

    do i = 1, size(self % ctrlPts)
      do j = 1, size(self % ctrlPts(i, :, :))
        if (self % ctrlPts(i, j, 1) < aabb(1)) aabb(1) = self % ctrlPts(i, j, 1) 
        if (self % ctrlPts(i, j, 2) < aabb(2)) aabb(2) = self % ctrlPts(i, j, 1) 
        if (self % ctrlPts(i, j, 1) > aabb(3)) aabb(3) = self % ctrlPts(i, j, 1) 
        if (self % ctrlPts(i, j, 2) > aabb(4)) aabb(4) = self % ctrlPts(i, j, 2) 
      end do 
    end do

    aabb(5) = self % minZ
    aabb(6) = self % maxZ


  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(bezierTwist), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    ! Parametric curve so no evaluation with global co-ords possible. 
    c = 0 

  end function evaluate


  !!
  !! Return distance to the surface
  !!
  !! Distance not available yet
  !!
  pure function distance(self, r, u) result(d)
    class(bezierTwist), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    ! Keep complier happy
    d = INF

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   For parallel direction halfspace is asigned by the sign of `evaluate` result.
  !!
  pure function going(self, r, u) result(halfspace)
    class(bezierTwist), intent(in)             :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace

    !! No going possible as surface tracking not yet supported. 
    halfspace = .false.

  end function going


  !!
  !! Checks if paticle is inside shape formed by bezier curves
  !!
  !! For c = F(r) halfspace is:
  !!   +1 : if outside cell
  !!   -1 : if inside cell
  !!
  !! Use surface tolerance if |c| = |F(r)| < SURF_TOL then halfspace is
  !! determined by direction of the particle (+ve if is is moving into +ve halfspace;
  !! -ve otherwise; determined with `going` procedure)
  !!
  !! Args:
  !!   r [in] -> Particle location
  !!   u [in] -> Particle direction (assume norm2(u) = 1.0 )
  !!
  !! Result:
  !!   True if position is in +ve halfspace. False if it is in -ve
  !!
  function halfspace(self, r, u) result(hs)
    class(bezierTwist), intent(in)                                    :: self
    real(defReal), dimension(3), intent(in)                           :: r
    real(defReal), dimension(3), intent(in)                           :: u
    real(defReal), dimension(3)                                       :: r_temp
    logical(defBool)                                                  :: hs
    real(defReal)                                                     :: t, cosine, sine, z, angle
    real(defReal), dimension(self % order + 1, 3)                     :: curve
    real(defReal), dimension(self % order + 1, 3)                     :: Q
    real(defReal), dimension(self % order + 1, 3)                     :: S
    real(defReal), dimension(3, 3)                                    :: rotationMatrix
    real(defReal)    :: x, y, zed
    real(defReal), dimension((self % order + 1) * (self % numCurves + 50), 3)       :: inner
    real(defReal), dimension((self % order + 1) * (self % numCurves + 50), 3)       :: outer
    real(defReal), dimension(self % numCurves + 50, self % order + 1, 3)  :: lastCtrlPts
    real(defReal), dimension(self % numCurves + 50, self % order + 1, 3)  :: newCtrlPts
    integer(shortInt)                                                 :: counter, i, j, m, n, a, b

    r_temp = r

    ! Keep compiler happy (in immpossible case of cell with no surfaces)
    hs = .false.

    ! if particle is out of upper and lower bounds then not in halfspace, return false
    if ((r(3) > self % maxZ) .or. (r(3) < self % minZ)) then
      hs = .false. 
      return
    end if

    ! Count number of subdivisions
    counter = 0
    n = 0
    m = 0
    a = 0 

    ! find rotation angle using z and twist period. 
    z= r(3)
    angle = modulo((z - self % minZ), self % twist) * 2.0 *  3.1415926535898 / self % twist
    
    ! Create rotation matrix depending on z value of particle and twist amount. 
    cosine = cos(angle)
    sine = sin(angle)

    
    ! Apply rotation to all initial control points
    do a = 1, self % numCurves
      do b = 1, self % order + 1
    
        x = self % ctrlPts(a, b, 1)
        y = self % ctrlPts(a, b, 2)

        ! Rotate points in the XY plane while keeping z-values unchanged
        lastCtrlPts(a, b, 1) = x * cosine - y * sine 
        lastCtrlPts(a, b, 2) = x * sine + y * cosine 
        lastCtrlPts(a, b, 3) = ZERO
        
      end do
    end do
    

    newCtrlPts = 0

    ! parametic subdivision point t
    t = 0.5

    main: do while (counter < 50)

      ! First calculate inner and outer paths within/around the bezier surfaces
      call self % findInnerOuterPaths(lastCtrlPts, inner, outer, counter, n, m)

      if (self % inOrOut(inner, n, r_temp)) then
        hs = .true.
        return
        
      elseif ( .not. self % inOrOut(outer, m, r_temp)) then
        hs = .false.
        return
      
      else

        ! if not inside either of the 
        do i = 1, self % numCurves + counter
          
          curve = lastCtrlPts(i, :, :)
          
          call self % order3Hull(curve)
  
          if (self % inConvexPoly(curve, r_temp)) then
            counter = counter + 1
            
            ! Subdivide curve at parametric point of t = 0.5
            call self % subdivide(curve, t, Q, S)
  
            ! Make new temp ctrl pts array - shuffle array along and make room for subdivided curve
            do j = 1, self % numCurves + counter
              if (j < i) then
                newCtrlPts(j, :, :) = lastCtrlPts(j, :, :)
              else if (j==i) then
                newCtrlPts(j, :, :) = Q(:, :)
              else if ( j == i + 1) then
                newCtrlPts(j, :, :) = S(:, :)
              else
                newCtrlPts(j, :, :) = lastCtrlPts(j-1, :, :)
              end if   
            end do
  
            lastCtrlPts = newCtrlPts
            
  
            ! go back to start of outside loop
            cycle main
          end if   
        end do
        
      
        ! Particle stuck and not found in any polygon
        ! Give particle nudge to set it on it's way
        r_temp = r_temp + (surf_tol * 10)
      
      
      end if 

    end do main




  end function halfspace


  !!
  !! Subdivision of Bezier Curve at point t
  !!
  !! Returns Bezier curves R and Q 
  !! 
  !! Currently Supports just cubic bezier curves, to be extended to higher orders
  !! Needs extending with loop to higher orders of bezier curve
  !!
  !! Uses De Casteljau’s Subdivision Algorithm with parametric 
  !! division point of t. Q is the LHS curve after split, and R is RHS curve
  !!
  subroutine subdivide(self, curve, t, Q, R)
    class(bezierTwist), intent(in)                               :: self
    real(defReal), dimension(self % order + 1, 3), intent(out)   :: Q
    real(defReal), dimension(self % order + 1, 3), intent(out)   :: R
    real(defReal), dimension(self % order + 1, 3), intent(in)    :: curve
    real(defReal), intent(in)                                    :: t
    real(defReal), dimension(3)                                  :: X  
    
    ! start and end points
    R(self % order + 1, :) = curve(self % order + 1, :)
    Q(1, :)                = curve(1, :) 

    ! X is the intermediattary point - not in either subdivided curve. 
    X    = (1-t) * curve(2, :) + t * curve(3, :) 
    R(3, :) = (1-t) * curve(3, :) + t * curve(4, :)
    Q(2, :) = (1-t) * curve(1, :) + t * curve(2, :)
    
    R(2, :) = (1-t) * X + t * R(3, :)
    Q(3, :) = (1-t) * Q(2, :) + t * X

    R(1, :) = (1-t) * Q(3, :) + t * R(2, :)
    Q(4, :) = (1-t) * Q(3, :) + t * R(2, :)

  end subroutine subdivide



  !!
  !! Find inner polygon and outer polygon around bezier curves
  !!
  !! Returns two lists of points which define inner and outer polygon
  !!
  !! Checks if point inside polygon and adjusts path acordingly
  !!
  subroutine findInnerOuterPaths(self, ctrlPts, inner, outer, counter, n, m)
    class(bezierTwist), intent(in)                                                  :: self
    real(defReal), dimension((self % order + 1) * (self % numCurves + 50), 3), intent(out)    :: inner
    real(defReal), dimension((self % order + 1) * (self % numCurves + 50), 3), intent(out)    :: outer
    real(defReal), dimension((self % order + 1) * (self % numCurves + 50), 3)                 :: poly
    real(defReal), dimension(self % numCurves + 50, self % order + 1, 3), intent(in):: ctrlPts
    integer(shortInt), intent(in)                                                   :: counter
    integer(shortInt), intent(out)                                                  :: m
    integer(shortInt), intent(out)                                                  :: n
    integer(shortInt)                                                               :: i, j

    poly = ctrlPts(:, 1, :)
    
    ! inner index = n
    ! outer index = m
    m = 1
    n = 1    
    do i = 1, self % numCurves + counter
      do j = 1, self % order
        !if (counter>1) print*, ctrlPts(i, j, :)
        if (j == 1) then
          inner(n, :) = ctrlPts(i, j, :)
          outer(m, :) = ctrlPts(i, j, :)
          !if (counter>1) print*, 'both'        
          n = n + 1
          m = m + 1
        else if (self % inOrOut(poly, self % numCurves + counter, ctrlPts(i, j, :))) then
          inner(n, :) = ctrlPts(i, j, :)
          !if (counter>1) print*, 'inner'
          n = n + 1
        else 
          outer(m, :) = ctrlPts(i, j, :)
          !if (counter>1) print*, 'outer'
          m = m + 1
        end if
      end do 
    end do 

    m = m - 1
    n = n - 1

  end subroutine findInnerOuterPaths

  


  !!
  !! In or out of aribitrary polygon. 
  !!
  !! Returns true or false wether in or out.   
  !!
  !! Doesnt matter if polygon is convex or concave, although if always convex,
  !! use other function as more efficient. Algorithm works via finding angles 
  !! between vectors from the point and summing in one direction. If angles 
  !! sum to zero then outside, if 2pi then point is inside.  
  !!
  !! Uses surface tolerance for if in or out, nudges if cant decide
  !!
  !! Currently uses 
  !!
  !!
  !!
  !!
  function inOrOut(self, path, len, r) result(inside)
    class(bezierTwist), intent(in)                                            :: self
    real(defReal), dimension((self % order + 1) * (self % numCurves + 50), 3), intent(in)  :: path
    real(defReal), dimension(3), intent(in)                                   :: r
    integer(shortInt), intent(in)                                             :: len
    real(defReal), dimension(2)                                               :: this_point, last_point
    real(defReal)                                                             :: total_angle, cosine_theta, theta
    real(defReal)                                                             :: cross, ONE
    logical(defBool)                                                          :: inside
    integer(shortInt)                                                         :: i

    inside = .false.

    ! index ranges added to ensure 2d aspect considered
    total_angle = 0
    ONE = 1

    last_point(1:2) = path(len, 1:2) - r(1:2) 

    do i = 1, len 
      this_point = path(i, 1:2) - r(1:2)
      if ((last_point(1) == 0) .and. (last_point(2) == 0)) print*, 'repeated points!', r(1:2), path(len, 1:2), len
      if ((this_point(1) == 0) .and. (this_point(2) == 0)) print*, 'repeated points!', r(1:2), path(i, 1:2), len

      cosine_theta = dot_product(last_point(1:2), this_point(1:2)) / ( norm2(this_point(1:2)) * norm2(last_point(1:2)) )
      
      ! resolve round off errors
      if (cosine_theta >=  1) cosine_theta =  0.99999999999
      if (cosine_theta <= -1) cosine_theta = -0.99999999999

      cross = last_point(1) * this_point(2) -  this_point(1) * last_point(2)

      theta = acos(cosine_theta) * sign(ONE, cross)
      
      total_angle = total_angle + theta

      last_point = this_point

      
    end do 


    if (abs(abs(total_angle) - 2 * 3.1415926535898) < surf_tol) then
      inside = .true.
      return

    else if (abs(total_angle) < surf_tol)  then
      inside = .false.
      return
    
    else
      print*, total_angle
      print*, r
      print *,'Bezier surface check failed'
    end if   

  end function inOrOut

  !!
  !! In or out of a convex polygon. 
  !!
  !! Returns true or false wether in or out.   
  !!
  !! Only works for convex polygons with number of points = order +1
  !! transforms polygon to be around the origin. 
  !!
  !!
  function inConvexPoly(self, poly, r) result(inside)
    class(bezierTwist), intent(in)                          :: self
    real(defReal), dimension(self % order + 1, 3), intent(in):: poly
    real(defReal), dimension(3), intent(in)                 :: r
    real(defReal), dimension(self % order + 1, 3)           :: poly_t
    real(defReal), dimension(2)                             :: last_point
    real(defReal)                                           :: last_angle, angle, ONE
    integer(shortInt)                                       :: i
    logical(defBool)                                        :: inside

    ! ONE = 1
    ONE = 1

    ! transform polygon by the point to the origin
    poly_t(:, 1) = poly(:, 1) - r(1)
    poly_t(:, 2) = poly(:, 2) - r(2)
    
    last_point = poly_t(self % order + 1, 1:2)
    last_angle = poly_t(1, 1) * last_point(2) - poly_t(1, 2) * last_point(1)  
    last_point = poly_t(1, 1:2)
    ! loop across angles to see if any |angle| > pi, which corresposnds to all "angles" being same sign 
    do i = 2, (self % order + 1)
      angle = poly_t(i, 1) * last_point(2) - poly_t(i, 2) * last_point(1)
      if (sign(ONE, angle) /= sign(ONE, last_angle)) then
        inside = .false.
        return
      end if
      last_point = poly_t(i, 1:2)
    end do

    ! If not returned false then point must be in polygon
    inside = .true.
    !print*, 'in convex hull of curves'

  end function inConvexPoly

  





  
  !!
  !! Simplified convex hull algorithm for bezier curves of order 3 
  !!
  !! Checks triangle combinations and turns directions  
  !! ie swap order of two points if crossed
  !! 
  !! To keep all arrays uniform, if the convex hull is a triangle, 
  !! we wrap the last point around again to fill the array
  !!
  subroutine order3Hull(self, curve)
    class(bezierTwist), intent(in)                  :: self
    real(defReal), intent(inout), dimension(4, 3)   :: curve
    real(defReal)                                   :: ABC, ABD, BCD, CAD, ONE
    real(defReal), dimension(3)                     :: A, B, C, D

    ONE=1

    
    A = curve(1, :)
    B = curve(2, :)
    C = curve(3, :)
    D = curve(4, :)

    ABC= sign(ONE, (A(2)-B(2))*C(1) + (B(1)-A(1))*C(2) + (A(1)*B(2)-B(1)*A(2)) )
    ABD= sign(ONE, (A(2)-B(2))*D(1) + (B(1)-A(1))*D(2) + (A(1)*B(2)-B(1)*A(2)) )
    BCD= sign(ONE, (B(2)-C(2))*D(1) + (C(1)-B(1))*D(2) + (B(1)*C(2)-C(1)*B(2)) )
    CAD= sign(ONE, (C(2)-A(2))*D(1) + (A(1)-C(1))*D(2) + (C(1)*A(2)-A(1)*C(2)) )

    if (ABC == ABD) then
      curve(1, :) = A
      curve(2, :) = B
      if (ABC == BCD) then 
        curve(3, :) = C
        if (ABC == CAD) then
          curve(4, :) = A
        else 
          curve(4, :) = D
        end if
      else
        curve(3, :) = D
        if (ABC == CAD) then
          curve(4, :) = C
        else
          curve(4, :) = A
        end if
      end if 
    else if (BCD == CAD) then
      curve(1, :) = A
      curve(2, :) = D
      curve(3, :) = B
      curve(4, :) = C
    else if (ABC == BCD) then
      curve(1, :) = B
      curve(2, :) = C
      curve(3, :) = D
      curve(4, :) = B
    else if (ABC == CAD) then
      curve(1, :) = C
      curve(2, :) = A
      curve(3, :) = D
      curve(4, :) = C
    else
      !print *,'3rd order convex hull not found' 
    end if

  end subroutine order3Hull




  !!
  !! Graham scan algorithm for finding convex hull of a bezier curve. 
  !!
  !! Currently not in use, see simplified 2d 3rd order algo above
  !!
  !!
  !! In future possible to add shortcut for curves of order 3, 
  !! ie swap order of two points if crossed
  !!
  subroutine grahamScan(self, curve)
    class(bezierTwist), intent(in)                              :: self
    real(defReal), intent(inout), dimension(self % order + 1, 3) :: curve


  end subroutine grahamScan




  !!
  !! Orientation of three point combination.
  !!
  !! Determines wether turn is Cw or CCW when finding convex hull
  !!
  !! This is then used during a Graham scan
  !!
  !! Currently implementation for 2D
  !!
  function orientation(self, p, q, r) result(ori)
    class(bezierTwist), intent(in)                           :: self
    real(defReal), intent(in), dimension(3)                 :: p
    real(defReal), intent(in), dimension(3)                 :: q
    real(defReal), intent(in), dimension(3)                 :: r
    real(defReal)                                           :: val
    Integer(shortInt)                                       :: ori

    val = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2))

    if (val == 0) then
      ori = 0
    elseif (val>0) then
      ori = 1
    else
      ori = -1
    end if

  end function orientation

  

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(bezierTwist), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % norm = ZERO
    self % offset = ZERO

  end subroutine kill

end module bezierTwist_class
