module trainedMLP_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  public :: trainedMLP
  public :: normaliseToBox
  public :: ACTIVATION_LEAKYRELU, ACTIVATION_RELU, ACTIVATION_TANH
  public :: MLP_MAX_DIM

  !! Activation type flags (stored in the weight file header)
  integer(shortInt), parameter :: ACTIVATION_LEAKYRELU = 1_shortInt
  integer(shortInt), parameter :: ACTIVATION_RELU      = 2_shortInt
  integer(shortInt), parameter :: ACTIVATION_TANH      = 3_shortInt

  !!
  !! Maximum layer width supported by the pure forward pass.
  !!
  !! evaluateRaw() uses stack-allocated work arrays of this size to stay pure
  !! (no ALLOCATE in a pure procedure). It bounds both the hidden width and the
  !! assembled input width. Increase if larger architectures are needed.
  !!
  integer(shortInt), parameter :: MLP_MAX_DIM = 512_shortInt

  !!
  !! Binary weight file magic number: ASCII "NSDF" = 0x4E534446, written and
  !! read as a 4-byte little-endian integer (native x86 byte order).
  !!
  integer(shortInt), parameter :: MAGIC_NUMBER   = int(z'4E534446', shortInt)
  integer(shortInt), parameter :: FORMAT_VERSION = 1_shortInt

  !! Tolerance for the embedded test-vector validation check
  real(defReal), parameter :: VALIDATION_TOL = 1.0e-10_defReal

  !!
  !! Trained MLP for signed distance function evaluation
  !!
  !! Holds architecture parameters and weight arrays for an MLP that maps
  !! normalised coordinates to a signed distance value. The forward pass
  !! (evaluate / evaluateRaw) is declared pure so it can be called from within
  !! SCONE's pure surface evaluate() procedure.
  !!
  !! This is shared machinery: neuralSurface_class holds one instance as its
  !! global SDF network, and deepLSSurface_class holds one instance as the
  !! shared DeepLS decoder (inputDim = latentDim + 3). Neither surface is a
  !! natural home for the type, hence the standalone module.
  !!
  !! Architecture (stored in the weight file header):
  !!   Layer 1:        inputDim  -> hiddenDim  (hidden activation)
  !!   Layers 2..N-1:  hiddenDim -> hiddenDim  (hidden activation)
  !!   Layer N:        hiddenDim -> 1          (tanh * sdfScale)
  !!
  !! Weight storage layout:
  !!   weights(i, j, l) = element (row i, col j) of the layer-l weight matrix.
  !!     Shape (maxDim, maxDim, numLayers) with maxDim = max(inputDim, hiddenDim).
  !!   biases(i, l)     = bias i for layer l.  Shape (maxDim, numLayers).
  !!
  !! Sign convention (evaluate / evaluateRaw): negative return = inside the
  !! surface (negative halfspace), positive return = outside.
  !!
  !! Public Members:
  !!   inputDim       -> Number of input dimensions fed to layer 1
  !!   hiddenDim      -> Width of every hidden layer
  !!   numLayers      -> Total number of weight matrices (>= 2)
  !!   activationType -> Hidden-layer activation (ACTIVATION_* flag)
  !!   leakyAlpha     -> LeakyReLU negative slope
  !!   sdfScale       -> tanh output multiplier
  !!   bboxMin        -> Minimum corner of the coordinate-normalisation box
  !!   bboxMax        -> Maximum corner of the coordinate-normalisation box
  !!   weights        -> Weight matrices; allocated by init(), filled by load()
  !!                     or written directly by the caller (see deepLSSurface)
  !!   biases         -> Bias vectors; same lifecycle as weights
  !!   isInit         -> True once init() has allocated the weight arrays
  !!
  !! Interface:
  !!   init        -> Allocate weight arrays for a given architecture
  !!   load        -> Read architecture + weights from a binary/text weight file
  !!   dump        -> Write architecture + weights to a binary weight file
  !!   evaluate    -> Pure forward pass on a world point (bbox-normalised here)
  !!   evaluateRaw -> Pure forward pass on a caller-assembled input vector
  !!   layerDims   -> Input/output dimensions of a given layer
  !!   kill        -> Deallocate weight arrays and reset to the uninitialised state
  !!
  type, public :: trainedMLP
    integer(shortInt) :: inputDim       = 0
    integer(shortInt) :: hiddenDim      = 0
    integer(shortInt) :: numLayers      = 0
    integer(shortInt) :: activationType = ACTIVATION_LEAKYRELU
    real(defReal)     :: leakyAlpha     = 0.01_defReal
    real(defReal)     :: sdfScale       = ONE

    !! Bounding box for coordinate normalisation: input coordinates are mapped
    !! onto [-1, 1] per axis before evaluation.
    real(defReal), dimension(3) :: bboxMin = ZERO
    real(defReal), dimension(3) :: bboxMax = ONE

    !! Weight matrices, allocated to the exact architecture dimensions at init()
    real(defReal), allocatable :: weights(:,:,:)
    real(defReal), allocatable :: biases(:,:)

    logical(defBool) :: isInit = .false.

  contains
    procedure :: init        => initMLP
    procedure :: load        => loadMLP
    procedure :: dump        => dumpMLP
    procedure :: evaluate    => evaluateMLP
    procedure :: evaluateRaw => evaluateRawMLP
    procedure :: layerDims   => layerDimsMLP
    procedure :: kill        => killMLP
  end type trainedMLP

contains

  !!
  !! Map a coordinate onto [-1, 1] given the normalisation box corners
  !!
  !! xn = 2 * (x - lo) / (hi - lo) - 1
  !!
  !! Elemental: works on scalars and, elementwise, on the dimension(3)
  !! coordinate and box-corner arrays used by the surfaces and the DeepLS
  !! reader alike.
  !!
  !! Args:
  !!   x  [in] -> coordinate value
  !!   lo [in] -> box minimum for this axis
  !!   hi [in] -> box maximum for this axis
  !!
  !! Result:
  !!   Normalised coordinate in [-1, 1] for x in [lo, hi].
  !!
  elemental function normaliseToBox(x, lo, hi) result(xn)
    real(defReal), intent(in) :: x, lo, hi
    real(defReal)             :: xn

    xn = TWO * (x - lo) / (hi - lo) - ONE

  end function normaliseToBox

  !!
  !! Allocate weight arrays for the given MLP architecture
  !!
  !! Must be called before weights are loaded. After this call the weights and
  !! biases arrays are allocated and zeroed.
  !!
  !! Args:
  !!   inputDim       [in] -> Number of input dimensions (3 for xyz)
  !!   hiddenDim      [in] -> Width of every hidden layer
  !!   numLayers      [in] -> Total number of weight matrices (>= 2)
  !!   activationType [in] -> Hidden-layer activation (ACTIVATION_* flag)
  !!   leakyAlpha     [in] -> LeakyReLU negative slope
  !!   sdfScale       [in] -> tanh output multiplier (maps [-1,1] to the SDF range)
  !!   bboxMin        [in] -> Minimum corner of the training bounding box
  !!   bboxMax        [in] -> Maximum corner of the training bounding box
  !!
  !! Errors:
  !!   fatalError if a dimension is invalid or hiddenDim > MLP_MAX_DIM
  !!
  subroutine initMLP(self, inputDim, hiddenDim, numLayers, activationType, &
                     leakyAlpha, sdfScale, bboxMin, bboxMax)
    class(trainedMLP), intent(inout)        :: self
    integer(shortInt), intent(in)           :: inputDim, hiddenDim, numLayers, activationType
    real(defReal), intent(in)               :: leakyAlpha, sdfScale
    real(defReal), dimension(3), intent(in) :: bboxMin, bboxMax
    integer(shortInt)                       :: maxDim
    character(100), parameter :: Here = 'initMLP (trainedMLP_class.f90)'

    ! Validate architecture parameters
    if (inputDim <= 0)           call fatalError(Here, 'inputDim must be positive')
    if (hiddenDim <= 0)          call fatalError(Here, 'hiddenDim must be positive')
    if (numLayers < 2)           call fatalError(Here, 'numLayers must be >= 2')
    if (hiddenDim > MLP_MAX_DIM) call fatalError(Here, 'hiddenDim exceeds MLP_MAX_DIM')
    if (sdfScale <= ZERO)        call fatalError(Here, 'sdfScale must be positive')

    ! Store architecture parameters
    self % inputDim       = inputDim
    self % hiddenDim      = hiddenDim
    self % numLayers      = numLayers
    self % activationType = activationType
    self % leakyAlpha     = leakyAlpha
    self % sdfScale       = sdfScale
    self % bboxMin        = bboxMin
    self % bboxMax        = bboxMax

    ! Allocate weight arrays padded to maxDim, which must accommodate both the
    ! input->hidden and hidden->hidden transitions.
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
  !! Input and output dimensions of layer l
  !!
  !! Layer 1 maps inputDim -> hiddenDim, interior layers hiddenDim -> hiddenDim,
  !! and the output layer hiddenDim -> 1.
  !!
  !! Args:
  !!   l      [in]  -> Layer index (1 .. numLayers)
  !!   inDim  [out] -> Number of inputs to layer l
  !!   outDim [out] -> Number of outputs from layer l
  !!
  pure subroutine layerDimsMLP(self, l, inDim, outDim)
    class(trainedMLP), intent(in)  :: self
    integer(shortInt), intent(in)  :: l
    integer(shortInt), intent(out) :: inDim, outDim

    if (l == 1_shortInt) then
      inDim  = self % inputDim
      outDim = self % hiddenDim
    else if (l < self % numLayers) then
      inDim  = self % hiddenDim
      outDim = self % hiddenDim
    else
      inDim  = self % hiddenDim
      outDim = 1_shortInt
    end if

  end subroutine layerDimsMLP

  !!
  !! Evaluate the MLP forward pass at a 3D world point
  !!
  !! Normalises the point onto [-1, 1] per axis using the training bounding box,
  !! then delegates to evaluateRaw. Declared pure so it can be called from
  !! within SCONE's pure surface evaluate() procedure.
  !!
  !! Sign convention: negative return = inside (negative halfspace),
  !!                  positive return = outside (positive halfspace).
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
    real(defReal), dimension(3)             :: xNorm

    xNorm = normaliseToBox(point, self % bboxMin, self % bboxMax)
    sdf   = self % evaluateRaw(xNorm)

  end function evaluateMLP

  !!
  !! Evaluate the MLP forward pass on a caller-assembled input vector, with NO
  !! bbox normalisation
  !!
  !! For DeepLS shared-decoder inference: one trainedMLP instance holds the
  !! shared decoder weights (inputDim = latentDim + 3) and the caller
  !! concatenates [per-voxel latent code, normalised local xyz] before calling
  !! this. Normalisation uses the voxel's own local bbox, which is stored on
  !! the surface rather than on this shared instance. See deepLSSurface_class.
  !!
  !! Declared pure so it can be called from within SCONE's pure surface
  !! evaluate() procedure. Uses fixed-size stack arrays (MLP_MAX_DIM) to avoid
  !! allocation in a pure context. Same sign convention as evaluate.
  !!
  !! Args:
  !!   input [in] -> Pre-normalised / pre-concatenated input, size >= inputDim
  !!
  !! Result:
  !!   Approximate signed distance value in world units.
  !!
  pure function evaluateRawMLP(self, input) result(sdf)
    class(trainedMLP), intent(in)           :: self
    real(defReal), dimension(:), intent(in) :: input
    real(defReal)                           :: sdf
    real(defReal), dimension(MLP_MAX_DIM)   :: h, hNext
    integer(shortInt)                       :: l, inDim, outDim
    logical(defBool)                        :: isHidden

    h(1:self % inputDim) = input(1:self % inputDim)

    do l = 1, self % numLayers
      isHidden = l < self % numLayers
      call self % layerDims(l, inDim, outDim)

      ! Linear transformation: hNext = W * h + b
      hNext(1:outDim) = matmul(self % weights(1:outDim, 1:inDim, l), h(1:inDim)) &
                      + self % biases(1:outDim, l)

      ! Hidden-layer activation (the output layer is left linear here; its
      ! tanh is applied once at the end)
      if (isHidden) then
        select case (self % activationType)
          case (ACTIVATION_LEAKYRELU)
            where (hNext(1:outDim) < ZERO) hNext(1:outDim) = self % leakyAlpha * hNext(1:outDim)
          case (ACTIVATION_RELU)
            where (hNext(1:outDim) < ZERO) hNext(1:outDim) = ZERO
          case (ACTIVATION_TANH)
            hNext(1:outDim) = tanh(hNext(1:outDim))
        end select
      end if

      h(1:outDim) = hNext(1:outDim)
    end do

    ! Output activation: tanh scaled by sdfScale
    sdf = tanh(h(1)) * self % sdfScale

  end function evaluateRawMLP

  !!
  !! Deallocate the weight arrays and return to the uninitialised state
  !!
  elemental subroutine killMLP(self)
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

  ! ---------------------------------------------------------------------------
  ! Weight file I/O
  ! ---------------------------------------------------------------------------

  !!
  !! Read a trained MLP from a binary or text weight file
  !!
  !! Detects the format via the magic number in the first 4 bytes and falls
  !! back to plain text if it is absent. Validates the loaded weights against
  !! the file's embedded test vector.
  !!
  !! Binary format (stream, little-endian):
  !!   Header:        magic(i4), version(i4), inputDim(i4), hiddenDim(i4),
  !!                  numLayers(i4), activationType(i4),
  !!                  leakyAlpha(r8), sdfScale(r8)
  !!   Normalisation: bboxMin(r8x3), bboxMax(r8x3)
  !!   Validation:    testInput(r8x3), testOutput(r8)
  !!   Per layer l:   weight(out,in)(r8, Fortran col-major), bias(out)(r8)
  !!
  !! Text format: the same values in the same order, one per line, with '#'
  !! comment lines skipped.
  !!
  !! Args:
  !!   filename [in] -> Path to the weight file
  !!
  !! Errors:
  !!   fatalError if the file cannot be opened, is the wrong version, is
  !!   truncated, or fails the embedded test-vector check.
  !!
  subroutine loadMLP(self, filename)
    class(trainedMLP), intent(out) :: self
    character(*), intent(in)       :: filename
    integer(shortInt)              :: unit, stat
    integer(shortInt)              :: magic
    real(defReal), dimension(3)    :: testInput
    real(defReal)                  :: testOutput
    character(100), parameter :: Here = 'loadMLP (trainedMLP_class.f90)'

    ! Open as an unformatted stream to probe the magic number
    open(newunit = unit,          &
         file    = filename,      &
         access  = 'stream',      &
         form    = 'unformatted', &
         status  = 'old',         &
         action  = 'read',        &
         iostat  = stat)
    if (stat /= 0) call fatalError(Here, 'Cannot open weight file: '//trim(filename))

    read(unit, pos=1, iostat=stat) magic

    if (stat == 0 .and. magic == MAGIC_NUMBER) then
      ! Binary format: continue reading from byte 5 onward
      call readBinary(self, unit, testInput, testOutput, Here)
      close(unit)
    else
      ! Text format fallback
      close(unit)
      open(newunit = unit,        &
           file    = filename,    &
           form    = 'formatted', &
           status  = 'old',       &
           action  = 'read',      &
           iostat  = stat)
      if (stat /= 0) call fatalError(Here, 'Cannot open weight file as text: '//trim(filename))
      call readText(self, unit, testInput, testOutput, Here)
      close(unit)
    end if

    call validateTestVector(self, testInput, testOutput, Here)

  end subroutine loadMLP

  !!
  !! Write a trained MLP to a binary weight file
  !!
  !! Used for testing and round-trip verification. The caller supplies a test
  !! point and the expected MLP output there; both are embedded in the file and
  !! checked on read.
  !!
  !! Args:
  !!   filename   [in] -> Path to write
  !!   testInput  [in] -> Test point embedded for validation on read
  !!   testOutput [in] -> Expected evaluate(testInput) embedded for validation
  !!
  !! Errors:
  !!   fatalError if the file cannot be opened for writing.
  !!
  subroutine dumpMLP(self, filename, testInput, testOutput)
    class(trainedMLP), intent(in)           :: self
    character(*), intent(in)                :: filename
    real(defReal), dimension(3), intent(in) :: testInput
    real(defReal), intent(in)               :: testOutput
    integer(shortInt)                       :: unit, stat, l
    integer(shortInt)                       :: inDim, outDim
    character(100), parameter :: Here = 'dumpMLP (trainedMLP_class.f90)'

    open(newunit = unit,          &
         file    = filename,      &
         access  = 'stream',      &
         form    = 'unformatted', &
         status  = 'replace',     &
         action  = 'write',       &
         iostat  = stat)
    if (stat /= 0) call fatalError(Here, 'Cannot open file for writing: '//trim(filename))

    ! Header
    write(unit) MAGIC_NUMBER
    write(unit) FORMAT_VERSION
    write(unit) self % inputDim
    write(unit) self % hiddenDim
    write(unit) self % numLayers
    write(unit) self % activationType
    write(unit) self % leakyAlpha
    write(unit) self % sdfScale
    ! Normalisation
    write(unit) self % bboxMin
    write(unit) self % bboxMax
    ! Test vector
    write(unit) testInput
    write(unit) testOutput

    ! Weights and biases, layer by layer
    do l = 1_shortInt, self % numLayers
      call self % layerDims(l, inDim, outDim)
      write(unit) self % weights(1:outDim, 1:inDim, l)
      write(unit) self % biases(1:outDim, l)
    end do

    close(unit)

  end subroutine dumpMLP

  ! ---------------------------------------------------------------------------
  ! Private I/O helpers
  ! ---------------------------------------------------------------------------

  !!
  !! Read the binary format from unit (position already at byte 5, past the magic)
  !!
  subroutine readBinary(mlp, unit, testInput, testOutput, Here)
    class(trainedMLP), intent(out)           :: mlp
    integer(shortInt), intent(in)            :: unit
    real(defReal), dimension(3), intent(out) :: testInput
    real(defReal), intent(out)               :: testOutput
    character(*), intent(in)                 :: Here
    integer(shortInt)           :: version, inputDim, hiddenDim, numLayers, activationType
    integer(shortInt)           :: l, inDim, outDim
    real(defReal)               :: leakyAlpha, sdfScale
    real(defReal), dimension(3) :: bboxMin, bboxMax
    integer(shortInt)           :: stat

    read(unit, pos=5, iostat=stat) version;        call ioCheck(stat, 'version', Here)
    read(unit,        iostat=stat) inputDim;       call ioCheck(stat, 'inputDim', Here)
    read(unit,        iostat=stat) hiddenDim;      call ioCheck(stat, 'hiddenDim', Here)
    read(unit,        iostat=stat) numLayers;      call ioCheck(stat, 'numLayers', Here)
    read(unit,        iostat=stat) activationType; call ioCheck(stat, 'activationType', Here)
    read(unit,        iostat=stat) leakyAlpha;     call ioCheck(stat, 'leakyAlpha', Here)
    read(unit,        iostat=stat) sdfScale;       call ioCheck(stat, 'sdfScale', Here)
    read(unit,        iostat=stat) bboxMin;        call ioCheck(stat, 'bboxMin', Here)
    read(unit,        iostat=stat) bboxMax;        call ioCheck(stat, 'bboxMax', Here)
    read(unit,        iostat=stat) testInput;      call ioCheck(stat, 'testInput', Here)
    read(unit,        iostat=stat) testOutput;     call ioCheck(stat, 'testOutput', Here)

    if (version /= FORMAT_VERSION) call fatalError(Here, 'Unsupported weight file version')

    call mlp % init(inputDim, hiddenDim, numLayers, activationType, &
                    leakyAlpha, sdfScale, bboxMin, bboxMax)

    do l = 1_shortInt, mlp % numLayers
      call mlp % layerDims(l, inDim, outDim)
      read(unit, iostat=stat) mlp % weights(1:outDim, 1:inDim, l)
      call ioCheck(stat, 'weights', Here)
      read(unit, iostat=stat) mlp % biases(1:outDim, l)
      call ioCheck(stat, 'biases', Here)
    end do

  end subroutine readBinary

  !!
  !! Read the text format from unit (already opened as formatted)
  !!
  !! Values are one per line; blank lines and lines starting with '#' are skipped.
  !!
  subroutine readText(mlp, unit, testInput, testOutput, Here)
    class(trainedMLP), intent(out)           :: mlp
    integer(shortInt), intent(in)            :: unit
    real(defReal), dimension(3), intent(out) :: testInput
    real(defReal), intent(out)               :: testOutput
    character(*), intent(in)                 :: Here
    integer(shortInt)           :: version, inputDim, hiddenDim, numLayers, activationType
    integer(shortInt)           :: magicText, l, inDim, outDim, i, j
    real(defReal)               :: leakyAlpha, sdfScale
    real(defReal), dimension(3) :: bboxMin, bboxMax

    call readTxtI(unit, magicText, Here)
    if (magicText /= MAGIC_NUMBER) &
      call fatalError(Here, 'Text weight file: expected magic number not found')
    call readTxtI(unit, version,        Here)
    call readTxtI(unit, inputDim,       Here)
    call readTxtI(unit, hiddenDim,      Here)
    call readTxtI(unit, numLayers,      Here)
    call readTxtI(unit, activationType, Here)
    call readTxtR(unit, leakyAlpha,     Here)
    call readTxtR(unit, sdfScale,       Here)
    call readTxtR(unit, bboxMin(1),     Here)
    call readTxtR(unit, bboxMin(2),     Here)
    call readTxtR(unit, bboxMin(3),     Here)
    call readTxtR(unit, bboxMax(1),     Here)
    call readTxtR(unit, bboxMax(2),     Here)
    call readTxtR(unit, bboxMax(3),     Here)
    call readTxtR(unit, testInput(1),   Here)
    call readTxtR(unit, testInput(2),   Here)
    call readTxtR(unit, testInput(3),   Here)
    call readTxtR(unit, testOutput,     Here)

    if (version /= FORMAT_VERSION) call fatalError(Here, 'Unsupported weight file version')

    call mlp % init(inputDim, hiddenDim, numLayers, activationType, &
                    leakyAlpha, sdfScale, bboxMin, bboxMax)

    ! Weights in Fortran column-major order (i varies fastest)
    do l = 1_shortInt, mlp % numLayers
      call mlp % layerDims(l, inDim, outDim)
      do j = 1_shortInt, inDim
        do i = 1_shortInt, outDim
          call readTxtR(unit, mlp % weights(i, j, l), Here)
        end do
      end do
      do i = 1_shortInt, outDim
        call readTxtR(unit, mlp % biases(i, l), Here)
      end do
    end do

  end subroutine readText

  !!
  !! Evaluate the MLP at testInput and check it matches the expected output
  !!
  subroutine validateTestVector(mlp, testInput, testOutput, Here)
    class(trainedMLP), intent(in)           :: mlp
    real(defReal), dimension(3), intent(in) :: testInput
    real(defReal), intent(in)               :: testOutput
    character(*), intent(in)                :: Here
    real(defReal)                           :: computed, err
    character(30)                           :: errStr

    computed = mlp % evaluate(testInput)
    err      = abs(computed - testOutput)

    if (err > VALIDATION_TOL) then
      write(errStr, '(ES12.4)') err
      call fatalError(Here, 'Weight file validation FAILED: embedded test vector mismatch. &
                            &|computed - expected| = '//trim(adjustl(errStr))// &
                            '. Likely cause: byte order or matrix transpose error.')
    end if

  end subroutine validateTestVector

  !!
  !! Abort with fatalError if a stream I/O status is non-zero, naming the field
  !!
  subroutine ioCheck(stat, field, Here)
    integer(shortInt), intent(in) :: stat
    character(*), intent(in)      :: field, Here

    if (stat /= 0) call fatalError(Here, 'I/O error reading field: '//trim(field))

  end subroutine ioCheck

  !!
  !! Read one integer from a formatted unit, skipping blank and '#' comment lines
  !!
  subroutine readTxtI(unit, val, Here)
    integer(shortInt), intent(in)  :: unit
    integer(shortInt), intent(out) :: val
    character(*), intent(in)       :: Here
    character(256)                 :: line
    integer(shortInt)             :: stat

    do
      read(unit, '(A)', iostat=stat) line
      if (stat /= 0) call fatalError(Here, 'Unexpected end of text weight file (reading integer)')
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=stat) val
      if (stat /= 0) call fatalError(Here, 'Cannot parse integer in text weight file: '//trim(line))
      return
    end do

  end subroutine readTxtI

  !!
  !! Read one real from a formatted unit, skipping blank and '#' comment lines
  !!
  subroutine readTxtR(unit, val, Here)
    integer(shortInt), intent(in) :: unit
    real(defReal), intent(out)    :: val
    character(*), intent(in)      :: Here
    character(256)                :: line
    integer(shortInt)            :: stat

    do
      read(unit, '(A)', iostat=stat) line
      if (stat /= 0) call fatalError(Here, 'Unexpected end of text weight file (reading real)')
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=stat) val
      if (stat /= 0) call fatalError(Here, 'Cannot parse real in text weight file: '//trim(line))
      return
    end do

  end subroutine readTxtR

end module trainedMLP_class
