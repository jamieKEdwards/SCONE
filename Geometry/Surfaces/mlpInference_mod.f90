module mlpInference_mod

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !! Activation type flags (stored in weight file header)
  integer(shortInt), parameter, public :: ACTIVATION_LEAKYRELU = 1_shortInt
  integer(shortInt), parameter, public :: ACTIVATION_RELU       = 2_shortInt
  integer(shortInt), parameter, public :: ACTIVATION_TANH       = 3_shortInt

  !!
  !! Maximum hidden layer width supported by the pure evaluate() forward pass.
  !! evaluate() uses stack-allocated work arrays of this size to remain pure
  !! (no ALLOCATE in pure procedures). Increase if larger architectures are needed.
  !!
  integer(shortInt), parameter, public :: MLP_MAX_DIM = 512_shortInt

  !!
  !! Trained MLP for signed distance function evaluation
  !!
  !! Holds architecture parameters and weight arrays for an MLP that maps
  !! 3D coordinates to a signed distance value. The forward pass (evaluate)
  !! is declared pure so it can be called from within SCONE's pure surface
  !! evaluate() procedure.
  !!
  !! Architecture (stored in weight file header):
  !!   - Layer 1:           inputDim  -> hiddenDim  (hidden activation)
  !!   - Layers 2..N-1:     hiddenDim -> hiddenDim  (hidden activation)
  !!   - Layer N (output):  hiddenDim -> 1          (tanh * sdfScale)
  !!
  !! Weight storage layout:
  !!   weights(i, j, l) = element (row i, col j) of weight matrix for layer l
  !!   Shape: (maxDim, maxDim, numLayers) where maxDim = max(inputDim, hiddenDim)
  !!   biases(i, l)     = bias i for layer l
  !!   Shape: (maxDim, numLayers)
  !!
  !! Members set by mlpWeightIO_mod after init():
  !!   weights, biases
  !!
  !! Interface:
  !!   init     -> Allocate weight arrays for given architecture
  !!   evaluate -> Pure forward pass; returns approximate SDF value
  !!   kill     -> Deallocate weight arrays
  !!
  type, public :: trainedMLP
    integer(shortInt) :: inputDim       = 0
    integer(shortInt) :: hiddenDim      = 0
    integer(shortInt) :: numLayers      = 0
    integer(shortInt) :: activationType = ACTIVATION_LEAKYRELU
    real(defReal)     :: leakyAlpha     = 0.01_defReal
    real(defReal)     :: sdfScale       = ONE

    !! Bounding box for coordinate normalisation.
    !! Input coordinates are mapped to [-1, 1] per axis before evaluation.
    real(defReal), dimension(3) :: bboxMin = ZERO
    real(defReal), dimension(3) :: bboxMax = ONE

    !! Weight matrices, allocated to exact architecture dimensions at init().
    !! Populated by mlpWeightIO_mod.
    real(defReal), allocatable :: weights(:,:,:)
    real(defReal), allocatable :: biases(:,:)

    logical(defBool) :: isInit = .false.

  contains
    procedure :: init     => initMLP
    procedure :: evaluate => evaluateMLP
    procedure :: kill     => killMLP
  end type trainedMLP

contains

  !!
  !! Allocate weight arrays for the given MLP architecture
  !!
  !! Must be called before weights are loaded. After this call, the
  !! weights and biases arrays are allocated and zeroed, ready for
  !! mlpWeightIO_mod to populate them.
  !!
  !! Args:
  !!   inputDim       [in] -> Number of input dimensions (3 for xyz)
  !!   hiddenDim      [in] -> Width of all hidden layers
  !!   numLayers      [in] -> Total number of weight matrices (>=2)
  !!   activationType [in] -> Activation for hidden layers (ACTIVATION_* flags)
  !!   leakyAlpha     [in] -> LeakyReLU negative slope
  !!   sdfScale       [in] -> tanh output multiplier (maps [-1,1] to SDF range)
  !!   bboxMin        [in] -> Minimum corner of training bounding box
  !!   bboxMax        [in] -> Maximum corner of training bounding box
  !!
  !! Errors:
  !!   fatalError if dimensions are invalid or hiddenDim > MLP_MAX_DIM
  !!
  subroutine initMLP(self, inputDim, hiddenDim, numLayers, activationType, &
                     leakyAlpha, sdfScale, bboxMin, bboxMax)
    class(trainedMLP), intent(inout)        :: self
    integer(shortInt), intent(in)           :: inputDim, hiddenDim, numLayers, activationType
    real(defReal), intent(in)               :: leakyAlpha, sdfScale
    real(defReal), dimension(3), intent(in) :: bboxMin, bboxMax
    integer(shortInt)                       :: maxDim
    character(100), parameter :: Here = 'initMLP (mlpInference_mod.f90)'

    ! Validate architecture parameters
    if (inputDim <= 0) then
      call fatalError(Here, 'inputDim must be positive')
    end if
    if (hiddenDim <= 0) then
      call fatalError(Here, 'hiddenDim must be positive')
    end if
    if (numLayers < 2) then
      call fatalError(Here, 'numLayers must be >= 2 (at least one hidden + output layer)')
    end if
    if (hiddenDim > MLP_MAX_DIM) then
      call fatalError(Here, 'hiddenDim exceeds MLP_MAX_DIM. Increase MLP_MAX_DIM to support this architecture.')
    end if
    if (sdfScale <= ZERO) then
      call fatalError(Here, 'sdfScale must be positive')
    end if

    ! Store architecture parameters
    self % inputDim       = inputDim
    self % hiddenDim      = hiddenDim
    self % numLayers      = numLayers
    self % activationType = activationType
    self % leakyAlpha     = leakyAlpha
    self % sdfScale       = sdfScale
    self % bboxMin        = bboxMin
    self % bboxMax        = bboxMax

    ! Allocate weight arrays padded to maxDim
    ! maxDim must accommodate both input->hidden and hidden->hidden transitions
    maxDim = max(inputDim, hiddenDim)
    if (allocated(self % weights)) deallocate(self % weights)
    if (allocated(self % biases))  deallocate(self % biases)
    allocate(self % weights(maxDim, maxDim, numLayers))
    allocate(self % biases(maxDim, numLayers))
    self % weights = ZERO
    self % biases  = ZERO

    self % isInit = .true.

  end subroutine initMLP

  !!
  !! Evaluate the MLP forward pass at a given 3D point
  !!
  !! Pure — safe to call from within SCONE's pure surface evaluate() procedure.
  !! Uses fixed-size stack arrays (MLP_MAX_DIM) to avoid allocation in pure context.
  !!
  !! Sign convention: negative return value = inside surface (negative halfspace),
  !!                  positive return value = outside surface (positive halfspace).
  !!
  !! Args:
  !!   point [in] -> 3D position in world coordinates
  !!
  !! Result:
  !!   Approximate signed distance value in world units.
  !!
  pure function evaluateMLP(self, point) result(sdf)
    class(trainedMLP), intent(in)           :: self
    real(defReal), dimension(3), intent(in) :: point
    real(defReal)                           :: sdf
    real(defReal), dimension(MLP_MAX_DIM)   :: h, h_next
    real(defReal), dimension(3)             :: x_norm
    integer(shortInt)                       :: l, i, inDim, outDim

    ! Normalise coordinates to [-1, 1] using training bounding box
    x_norm = TWO * (point - self % bboxMin) / (self % bboxMax - self % bboxMin) - ONE

    ! Load normalised input into work vector
    h(1:self % inputDim) = x_norm(1:self % inputDim)
    inDim = self % inputDim

    ! Forward pass through each layer
    do l = 1, self % numLayers

      ! Output dimension: hiddenDim for hidden layers, 1 for output layer
      if (l < self % numLayers) then
        outDim = self % hiddenDim
      else
        outDim = 1
      end if

      ! Linear transformation: h_next = W * h + b
      h_next(1:outDim) = matmul(self % weights(1:outDim, 1:inDim, l), h(1:inDim)) &
                       + self % biases(1:outDim, l)

      ! Apply hidden layer activation (not applied to output layer)
      if (l < self % numLayers) then
        select case (self % activationType)
          case (ACTIVATION_LEAKYRELU)
            do i = 1, outDim
              if (h_next(i) < ZERO) h_next(i) = self % leakyAlpha * h_next(i)
            end do
          case (ACTIVATION_RELU)
            do i = 1, outDim
              if (h_next(i) < ZERO) h_next(i) = ZERO
            end do
          case (ACTIVATION_TANH)
            h_next(1:outDim) = tanh(h_next(1:outDim))
        end select
      end if

      ! Advance to next layer
      h(1:outDim) = h_next(1:outDim)
      inDim = outDim

    end do

    ! Output activation: tanh scaled by sdfScale
    sdf = tanh(h(1)) * self % sdfScale

  end function evaluateMLP

  !!
  !! Return to uninitialised state and deallocate weight arrays
  !!
  subroutine killMLP(self)
    class(trainedMLP), intent(inout) :: self

    if (allocated(self % weights)) deallocate(self % weights)
    if (allocated(self % biases))  deallocate(self % biases)

    self % inputDim       = 0
    self % hiddenDim      = 0
    self % numLayers      = 0
    self % activationType = ACTIVATION_LEAKYRELU
    self % leakyAlpha     = 0.01_defReal
    self % sdfScale       = ONE
    self % bboxMin        = ZERO
    self % bboxMax        = ONE
    self % isInit         = .false.

  end subroutine killMLP

end module mlpInference_mod
