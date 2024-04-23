module NN
   use ftorch

   implicit none

   integer, parameter :: nsize = 10
   real(kind=8), save, dimension(nsize, nsize) :: fc1w
   real(kind=8), save, dimension(nsize, 1) :: fc1b
   real(kind=8), save, dimension(2, nsize) :: fc2w
   real(kind=8), save, dimension(2, 1) :: fc2b

contains

   subroutine read_wb(nrows, ncols, filename, lt)
      implicit none
      integer :: file_unit, status
      integer :: nrows, ncols, i, j
      character(len=100) :: filename
      character(len=4) :: lt

      ! filename = "./models/nn-02/w-b/fc1w.txt"
      ! nrows = 8
      ! ncols = 20
      open (unit=file_unit, file=trim(filename), status='old', action='read', iostat=status)
      if (status /= 0) then
         print *, "Error opening file ", trim(filename)
         stop
      end if

      print *, "LT:"
      print *, trim(lt)
      print *, "\n"
      if (trim(lt) == "fc1w") then
         do i = 1, nrows
            read (file_unit, *) (fc1w(i, j), j=1, ncols)
         end do
      else if (trim(lt) == "fc1b") then
         do i = 1, nrows
            read (file_unit, *) (fc1b(i, j), j=1, ncols)
         end do
      else if (trim(lt) == "fc2w") then
         do i = 1, nrows
            read (file_unit, *) (fc2w(i, j), j=1, ncols)
         end do
      else if (trim(lt) == "fc2b") then
         do i = 1, nrows
            read (file_unit, *) (fc2b(i, j), j=1, ncols)
         end do
      end if

      close (file_unit)

   end subroutine

   subroutine read_wb_reduced(nrows, ncols, filename, lt)
      implicit none
      integer :: file_unit, status
      integer :: nrows, ncols, i, j
      character(len=100) :: filename
      character(len=4) :: lt

      ! filename = "./models/nn-02/w-b/fc1w.txt"
      ! nrows = 8
      ! ncols = 20
      open (unit=file_unit, file=trim(filename), status='old', action='read', iostat=status)
      if (status /= 0) then
         print *, "Error opening file ", trim(filename)
         stop
      end if

      print *, "LT:"
      print *, trim(lt)
      print *, "\n"
      if (trim(lt) == "fc1w") then
         do i = 1, nrows
            read (file_unit, *) (fc1w(i, j), j=1, ncols)
         end do
      else if (trim(lt) == "fc1b") then
         do i = 1, nrows
            read (file_unit, *) (fc1b(i, j), j=1, ncols)
         end do
      else if (trim(lt) == "fc2w") then
         do i = 1, nrows
            read (file_unit, *) (fc2w(i, j), j=1, ncols)
         end do
      else if (trim(lt) == "fc2b") then
         do i = 1, nrows
            read (file_unit, *) (fc2b(i, j), j=1, ncols)
         end do
      end if

      close (file_unit)

   end subroutine


   subroutine l_relu(arr, N)
      implicit none

      integer, intent(in) :: N
      real(kind=8), intent(INOUT), dimension(N, 1) :: arr
      integer :: i

      do i = 1, N
         if (arr(i, 1) < 0.0) then
            arr(i, 1) = -0.01*arr(i, 1)
         end if
      end do

   end subroutine l_relu

   subroutine read_w(file_name, w_name, mat, N1, N2)
      real(kind=8), intent(out), dimension(:, :) :: mat
      character(len=*) :: file_name, w_name
      integer, intent(in) :: N1, N2

      open (unit=10, file=file_name, status="old")

      close (10)

   end subroutine read_w

   subroutine classify(input, classification)
      real(kind=8), intent(inout), dimension(nsize, 1) :: input
      integer, intent(out) :: classification

      ! real(kind=8), dimension(size, 20) :: fc1_w
      ! real(kind=8), dimension(20, 1) :: fc1_b

      ! real(kind=8), dimension(2, 20) :: fc2_w
      ! real(kind=8), dimension(2, 1) :: fc2_b

      real(kind=8), dimension(nsize, 1) :: temp1
      real(kind=8), dimension(2, 1) :: temp2

      real(kind=8) :: max_val

      temp1 = fc1b + matmul(fc1w, input)
      call l_relu(temp1, nsize)
      temp2 = fc2b + matmul(fc2w, temp1)

      temp2(:, 1) = temp2(:, 1) - maxval(temp2(:, 1))
      if (exp(temp2(2, 1))/(exp(temp2(1, 1)) + exp(temp2(2, 1))) >= 0.99) then
         classification = 1
      else
         classification = 0
      end if

      !F.softmax(fc2_b.unsqueeze(-1) + (fc2_w @ F.leaky_relu(fc1_b.unsqueeze(-1) + fc1_w @ test.unsqueeze(-1))),dim=0)

   end subroutine classify

   subroutine classify_ftorch(input_in, classification)
      real(kind=4), intent(inout), dimension(20, 1) :: input_in
      integer, intent(out) :: classification

      type(torch_module) :: model
      integer, parameter :: n_inputs = 1
      type(torch_tensor), dimension(n_inputs) :: model_input_arr
      type(torch_tensor) :: model_output
      real(kind=4), dimension(1, 20), target  :: input
      real(kind=4), dimension(1, 2), target   :: output

      ! Set up number of dimensions of input tensor and axis order
      integer, parameter :: in_dims = 2
      integer :: in_layout(in_dims) = [1, 2]
      integer, parameter :: out_dims = 2
      integer :: out_layout(out_dims) = [1, 2]

      model = torch_module_load("./models/nn-02/script-nn-02.pt")

      input(1, :) = input_in(:, 1)

      model_input_arr(1) = torch_tensor_from_array(input, in_layout, torch_kCPU)
      model_output = torch_tensor_from_array(output, out_layout, torch_kCPU)

      call torch_module_forward(model, model_input_arr, n_inputs, model_output)

      if (output(1, 1) >= output(1, 2)) then
         classification = 0
      else
         classification = 1
      end if

      call torch_module_delete(model)
      call torch_tensor_delete(model_input_arr(1))
      call torch_tensor_delete(model_output)

   end subroutine classify_ftorch

end module NN
