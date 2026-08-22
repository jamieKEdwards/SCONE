module deepLSSurface_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface
  use mlpInference_mod,   only : trainedMLP
  use deepLSWeightIO_mod, only : readDeepLSWeights, &
                                  DEEPLS_CONST_INSIDE, DEEPLS_CONST_OUTSIDE, &
                                  DEEPLS_HAS_MLP

  implicit none
  private

  character(*), parameter :: TYPE_NAME = 'deepLSSurface'

  !!
  !! DeepLS surface, paper-faithful architecture (Chabra et al., ECCV 2020)
  !!
  !! ONE shared decoder network (same size as the global MLP: 4 layers, 128
  !! hidden by default) conditioned on a per-voxel latent code, rather than
  !! one global network representing the whole domain. The decoder is
  !! evaluated identically for every voxel; only the latent code and the
  !! local normalisation bbox change. Trained via
  !! scripts/neuralSurface/train_deepls.py.
  !!
  !! Sign convention identical to neuralSurface: F(r) < 0 => inside,
  !! F(r) > 0 => outside. Delta (Woodcock) tracking only.
  !!
  !! Private Members:
  !!   decoder      -> the ONE shared trainedMLP (inputDim = latentDim + 3)
  !!   latentDim    -> latent code length
  !!   nvox/gridOrigin/voxelSize -> voxel grid geometry
  !!   voxelStatus  -> (nx,ny,nz) DEEPLS_* status codes
  !!   voxelBboxMin/Max -> (3,nx,ny,nz) per-voxel local normalisation bbox
  !!   voxelLatent  -> (latentDim,nx,ny,nz) per-voxel latent code
  !!   geomScale    -> divides world-space coordinates before the voxel-grid
  !!                  lookup, letting one trained weight file represent
  !!                  geometrically similar surfaces of a different physical
  !!                  size without retraining
  !!
  !! Sample Dictionary Input:
  !!   surf { type deepLSSurface; id 1; weightFile ./weights.bin;
  !!          # geometricScale 1.0; #
  !!        }
  !!
  !! See misclassClerk_class (Tallies/TallyClerks) for a halfspace
  !! misclassification diagnostic against a reference region.
  !!
  type, public, extends(surface) :: deepLSSurface
    private
    type(trainedMLP)                                 :: decoder
    integer(shortInt)                                :: latentDim  = 0
    integer(shortInt), dimension(3)                  :: nvox       = 0
    real(defReal), dimension(3)                      :: gridOrigin = ZERO
    real(defReal), dimension(3)                      :: voxelSize  = ONE
    integer(shortInt), dimension(:,:,:), allocatable :: voxelStatus
    real(defReal), dimension(:,:,:,:), allocatable   :: voxelBboxMin, voxelBboxMax
    real(defReal), dimension(:,:,:,:), allocatable   :: voxelLatent
    real(defReal)                                    :: geomScale  = ONE
  contains
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
  end type deepLSSurface

contains

  pure function myType(self) result(str)
    class(deepLSSurface), intent(in) :: self
    character(:), allocatable        :: str

    str = TYPE_NAME

  end function myType

  !!
  !! Initialise from dictionary
  !!
  !! See surface_inter for details
  !!
  subroutine init(self, dict)
    class(deepLSSurface), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: id
    character(pathLen)                  :: weightFile
    character(100), parameter :: Here = 'init (deepLSSurface_class.f90)'

    call dict % get(id, 'id')
    if (id < 1) call fatalError(Here, 'Invalid surface id: '//numToChar(id))

    call dict % get(weightFile, 'weightFile')

    call readDeepLSWeights(self % decoder, self % latentDim, self % nvox, &
                           self % gridOrigin, self % voxelSize, self % voxelStatus, &
                           self % voxelBboxMin, self % voxelBboxMax, self % voxelLatent, &
                           trim(weightFile))
    call self % setId(id)

    if (dict % isPresent('geometricScale')) then
      call dict % get(self % geomScale, 'geometricScale')
      if (self % geomScale <= ZERO) &
        call fatalError(Here, 'geometricScale must be positive')
    else
      self % geomScale = ONE
    end if

  end subroutine init

  pure function boundingBox(self) result(aabb)
    class(deepLSSurface), intent(in) :: self
    real(defReal), dimension(6)      :: aabb

    aabb(1:3) = self % gridOrigin * self % geomScale
    aabb(4:6) = (self % gridOrigin + self % nvox * self % voxelSize) * self % geomScale

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(deepLSSurface), intent(in)  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: rs, xn
    integer(shortInt), dimension(3)         :: idx
    real(defReal), dimension(self % latentDim + 3) :: decoderInput

    rs = r / self % geomScale

    if (any(rs < self % gridOrigin) .or. &
        any(rs > self % gridOrigin + self % nvox * self % voxelSize)) then
      c = ONE * self % geomScale
      return
    end if

    idx = floor((rs - self % gridOrigin) / self % voxelSize) + 1
    idx = max(1, min(self % nvox, idx))

    select case (self % voxelStatus(idx(1), idx(2), idx(3)))
      case (DEEPLS_CONST_INSIDE)
        c = -ONE
      case (DEEPLS_CONST_OUTSIDE)
        c = ONE
      case (DEEPLS_HAS_MLP)
        decoderInput(1:self % latentDim) = self % voxelLatent(:, idx(1), idx(2), idx(3))
        xn = TWO * (rs - self % voxelBboxMin(:, idx(1), idx(2), idx(3))) / &
             (self % voxelBboxMax(:, idx(1), idx(2), idx(3)) - &
              self % voxelBboxMin(:, idx(1), idx(2), idx(3))) - ONE
        decoderInput(self % latentDim + 1 : self % latentDim + 3) = xn
        c = self % decoder % evaluateRaw(decoderInput)
      case default
        c = ONE
    end select

    c = c * self % geomScale

  end function evaluate

  pure function distance(self, r, u) result(d)
    class(deepLSSurface), intent(in)  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d

    d = INF

  end function distance

  pure function going(self, r, u) result(hs)
    class(deepLSSurface), intent(in)  :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    real(defReal), parameter :: FD_STEP = 1.0e-7_defReal

    hs = self % evaluate(r + FD_STEP * u) > ZERO

  end function going

end module deepLSSurface_class
