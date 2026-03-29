module mlpWeightIO_mod

  use numPrecision
  use genericProcedures, only : fatalError
  use mlpInference_mod,  only : trainedMLP

  implicit none
  private

  !!
  !! Binary file magic number: ASCII "NSDF" = 0x4E534446
  !! Written and read as a 4-byte little-endian integer (native x86 byte order).
  !!
  integer(shortInt), parameter :: MAGIC_NUMBER   = int(z'4E534446', shortInt)
  integer(shortInt), parameter :: FORMAT_VERSION = 1_shortInt

  !!
  !! Tolerance for the embedded test vector validation check.
  !!
  real(defReal), parameter :: VALIDATION_TOL = 1.0e-10_defReal

  public :: readMLPWeights
  public :: writeMLPWeights

contains

  !!
  !! Read a trained MLP from a binary or text weight file
  !!
  !! Detects format via magic number in the first 4 bytes.
  !! Falls back to plain-text if magic number is absent.
  !! Validates weights after loading using the embedded test vector.
  !!
  !! Binary format (stream, little-endian):
  !!   Header:       magic(i4), version(i4), input_dim(i4), hidden_dim(i4),
  !!                 num_layers(i4), activation_type(i4),
  !!                 leaky_alpha(r8), sdf_scale(r8)
  !!   Normalisation: bbox_min(r8×3), bbox_max(r8×3)
  !!   Validation:   test_input(r8×3), test_output(r8)
  !!   Per layer l:  weight(out,in)(r8, Fortran col-major), bias(out)(r8)
  !!
  !! Text format: same values in same order, one per line, '#' comment lines skipped.
  !!
  !! Args:
  !!   mlp      [out] -> Loaded and validated trainedMLP
  !!   filename [in]  -> Path to weight file
  !!
  subroutine readMLPWeights(mlp, filename)
    type(trainedMLP), intent(out) :: mlp
    character(*), intent(in)      :: filename
    integer                       :: unit, stat
    integer(shortInt)             :: magic
    real(defReal), dimension(3)   :: test_input
    real(defReal)                 :: test_output
    character(100), parameter :: Here = 'readMLPWeights (mlpWeightIO_mod.f90)'

    ! Open as unformatted stream to probe magic number
    open(newunit = unit,           &
         file    = filename,       &
         access  = 'stream',       &
         form    = 'unformatted',  &
         status  = 'old',          &
         action  = 'read',         &
         iostat  = stat)
    if (stat /= 0) call fatalError(Here, 'Cannot open weight file: '//trim(filename))

    ! Read first 4 bytes
    read(unit, pos=1, iostat=stat) magic

    if (stat == 0 .and. magic == MAGIC_NUMBER) then
      ! Binary format: continue reading from byte 5 onward
      call readBinary(mlp, unit, test_input, test_output, Here)
      close(unit)
    else
      ! Text format fallback
      close(unit)
      open(newunit = unit,          &
           file    = filename,      &
           form    = 'formatted',   &
           status  = 'old',         &
           action  = 'read',        &
           iostat  = stat)
      if (stat /= 0) call fatalError(Here, 'Cannot open weight file as text: '//trim(filename))
      call readText(mlp, unit, test_input, test_output, Here)
      close(unit)
    end if

    ! Validate with embedded test vector
    call validateTestVector(mlp, test_input, test_output, Here)

  end subroutine readMLPWeights


  !!
  !! Write a trained MLP to a binary weight file
  !!
  !! Used for testing and for round-trip verification.
  !! The caller supplies a test point and the expected MLP output at that point,
  !! which are embedded in the file and checked on read.
  !!
  subroutine writeMLPWeights(mlp, filename, test_input, test_output)
    type(trainedMLP), intent(in)            :: mlp
    character(*), intent(in)               :: filename
    real(defReal), dimension(3), intent(in) :: test_input
    real(defReal), intent(in)              :: test_output
    integer                                :: unit, stat, l
    integer(shortInt)                      :: in_dim, out_dim
    character(100), parameter :: Here = 'writeMLPWeights (mlpWeightIO_mod.f90)'

    open(newunit = unit,           &
         file    = filename,       &
         access  = 'stream',       &
         form    = 'unformatted',  &
         status  = 'replace',      &
         action  = 'write',        &
         iostat  = stat)
    if (stat /= 0) call fatalError(Here, 'Cannot open file for writing: '//trim(filename))

    ! Header
    write(unit) MAGIC_NUMBER
    write(unit) FORMAT_VERSION
    write(unit) mlp % inputDim
    write(unit) mlp % hiddenDim
    write(unit) mlp % numLayers
    write(unit) mlp % activationType
    write(unit) mlp % leakyAlpha
    write(unit) mlp % sdfScale
    ! Normalisation
    write(unit) mlp % bboxMin
    write(unit) mlp % bboxMax
    ! Test vector
    write(unit) test_input
    write(unit) test_output

    ! Weights and biases, layer by layer
    do l = 1, mlp % numLayers
      call layerDims(mlp, int(l, shortInt), in_dim, out_dim)
      write(unit) mlp % weights(1:out_dim, 1:in_dim, l)
      write(unit) mlp % biases(1:out_dim, l)
    end do

    close(unit)

  end subroutine writeMLPWeights


  ! ---------------------------------------------------------------------------
  ! Private helpers
  ! ---------------------------------------------------------------------------

  !!
  !! Read binary format from unit (file position already at byte 5 after magic read)
  !!
  subroutine readBinary(mlp, unit, test_input, test_output, Here)
    type(trainedMLP), intent(out)            :: mlp
    integer, intent(in)                      :: unit
    real(defReal), dimension(3), intent(out) :: test_input
    real(defReal), intent(out)               :: test_output
    character(*), intent(in)                 :: Here
    integer(shortInt) :: version, input_dim, hidden_dim, num_layers, activation_type
    integer(shortInt) :: l, in_dim, out_dim
    real(defReal)     :: leaky_alpha, sdf_scale
    real(defReal), dimension(3) :: bbox_min, bbox_max
    integer           :: stat

    ! Read from byte 5 (pos=5 is right after the 4-byte magic number)
    read(unit, pos=5,  iostat=stat) version;          call ioCheck(stat, 'version', Here)
    read(unit,         iostat=stat) input_dim;         call ioCheck(stat, 'input_dim', Here)
    read(unit,         iostat=stat) hidden_dim;        call ioCheck(stat, 'hidden_dim', Here)
    read(unit,         iostat=stat) num_layers;        call ioCheck(stat, 'num_layers', Here)
    read(unit,         iostat=stat) activation_type;   call ioCheck(stat, 'activation_type', Here)
    read(unit,         iostat=stat) leaky_alpha;       call ioCheck(stat, 'leaky_alpha', Here)
    read(unit,         iostat=stat) sdf_scale;         call ioCheck(stat, 'sdf_scale', Here)
    read(unit,         iostat=stat) bbox_min;          call ioCheck(stat, 'bbox_min', Here)
    read(unit,         iostat=stat) bbox_max;          call ioCheck(stat, 'bbox_max', Here)
    read(unit,         iostat=stat) test_input;        call ioCheck(stat, 'test_input', Here)
    read(unit,         iostat=stat) test_output;       call ioCheck(stat, 'test_output', Here)

    if (version /= FORMAT_VERSION) call fatalError(Here, 'Unsupported weight file version')

    call mlp % init(input_dim, hidden_dim, num_layers, activation_type, &
                    leaky_alpha, sdf_scale, bbox_min, bbox_max)

    do l = 1_shortInt, mlp % numLayers
      call layerDims(mlp, l, in_dim, out_dim)
      read(unit, iostat=stat) mlp % weights(1:out_dim, 1:in_dim, l)
      call ioCheck(stat, 'weights', Here)
      read(unit, iostat=stat) mlp % biases(1:out_dim, l)
      call ioCheck(stat, 'biases', Here)
    end do

  end subroutine readBinary


  !!
  !! Read text format from unit (already opened as formatted)
  !! Values are one per line; lines starting with '#' are skipped.
  !!
  subroutine readText(mlp, unit, test_input, test_output, Here)
    type(trainedMLP), intent(out)            :: mlp
    integer, intent(in)                      :: unit
    real(defReal), dimension(3), intent(out) :: test_input
    real(defReal), intent(out)               :: test_output
    character(*), intent(in)                 :: Here
    integer(shortInt) :: version, input_dim, hidden_dim, num_layers, activation_type
    integer(shortInt) :: magic_text, l, in_dim, out_dim
    integer(shortInt) :: i, j
    real(defReal)     :: leaky_alpha, sdf_scale
    real(defReal), dimension(3) :: bbox_min, bbox_max

    call readTxtI(unit, magic_text,       Here)
    if (magic_text /= MAGIC_NUMBER) &
      call fatalError(Here, 'Text weight file: expected magic number not found')
    call readTxtI(unit, version,          Here)
    call readTxtI(unit, input_dim,        Here)
    call readTxtI(unit, hidden_dim,       Here)
    call readTxtI(unit, num_layers,       Here)
    call readTxtI(unit, activation_type,  Here)
    call readTxtR(unit, leaky_alpha,      Here)
    call readTxtR(unit, sdf_scale,        Here)
    call readTxtR(unit, bbox_min(1),      Here)
    call readTxtR(unit, bbox_min(2),      Here)
    call readTxtR(unit, bbox_min(3),      Here)
    call readTxtR(unit, bbox_max(1),      Here)
    call readTxtR(unit, bbox_max(2),      Here)
    call readTxtR(unit, bbox_max(3),      Here)
    call readTxtR(unit, test_input(1),    Here)
    call readTxtR(unit, test_input(2),    Here)
    call readTxtR(unit, test_input(3),    Here)
    call readTxtR(unit, test_output,      Here)

    if (version /= FORMAT_VERSION) call fatalError(Here, 'Unsupported weight file version')

    call mlp % init(input_dim, hidden_dim, num_layers, activation_type, &
                    leaky_alpha, sdf_scale, bbox_min, bbox_max)

    ! Weights in Fortran column-major order (i varies fastest)
    do l = 1_shortInt, mlp % numLayers
      call layerDims(mlp, l, in_dim, out_dim)
      do j = 1_shortInt, in_dim
        do i = 1_shortInt, out_dim
          call readTxtR(unit, mlp % weights(i, j, l), Here)
        end do
      end do
      do i = 1_shortInt, out_dim
        call readTxtR(unit, mlp % biases(i, l), Here)
      end do
    end do

  end subroutine readText


  !!
  !! Evaluate the MLP at test_input and check it matches expected output
  !!
  subroutine validateTestVector(mlp, test_input, test_output, Here)
    type(trainedMLP), intent(in)            :: mlp
    real(defReal), dimension(3), intent(in) :: test_input
    real(defReal), intent(in)               :: test_output
    character(*), intent(in)                :: Here
    real(defReal)                           :: computed, err
    character(30)                           :: errStr

    computed = mlp % evaluate(test_input)
    err      = abs(computed - test_output)

    if (err > VALIDATION_TOL) then
      write(errStr, '(ES12.4)') err
      call fatalError(Here, 'Weight file validation FAILED: embedded test vector mismatch. &
                            &|computed - expected| = '//trim(adjustl(errStr))// &
                            '. Likely cause: byte order or matrix transpose error.')
    end if

  end subroutine validateTestVector


  !!
  !! Return input and output dimensions for layer l
  !!
  pure subroutine layerDims(mlp, l, in_dim, out_dim)
    type(trainedMLP), intent(in)   :: mlp
    integer(shortInt), intent(in)  :: l
    integer(shortInt), intent(out) :: in_dim, out_dim

    if (l == 1_shortInt) then
      in_dim  = mlp % inputDim
      out_dim = mlp % hiddenDim
    else if (l < mlp % numLayers) then
      in_dim  = mlp % hiddenDim
      out_dim = mlp % hiddenDim
    else
      in_dim  = mlp % hiddenDim
      out_dim = 1_shortInt
    end if

  end subroutine layerDims


  !!
  !! Abort with fatalError if stat /= 0
  !!
  subroutine ioCheck(stat, field, Here)
    integer, intent(in)      :: stat
    character(*), intent(in) :: field, Here
    if (stat /= 0) call fatalError(Here, 'I/O error reading field: '//trim(field))
  end subroutine ioCheck


  !!
  !! Read one integer(shortInt) from a formatted unit, skipping blank lines and '#' comments
  !!
  subroutine readTxtI(unit, val, Here)
    integer, intent(in)            :: unit
    integer(shortInt), intent(out) :: val
    character(*), intent(in)       :: Here
    character(256)                 :: line
    integer                        :: stat
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
  !! Read one real(defReal) from a formatted unit, skipping blank lines and '#' comments
  !!
  subroutine readTxtR(unit, val, Here)
    integer, intent(in)        :: unit
    real(defReal), intent(out) :: val
    character(*), intent(in)   :: Here
    character(256)             :: line
    integer                    :: stat
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

end module mlpWeightIO_mod
