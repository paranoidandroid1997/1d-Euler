module NN
    use ftorch

    implicit none

    real(kind = 8), save, dimension(20,20) :: fc1w
    real(kind = 8), save, dimension(20,1) :: fc1b
    real(kind = 8), save, dimension(2,20) :: fc2w
    real(kind = 8), save, dimension(2,1) :: fc2b


contains

    subroutine read_wb(nrows, ncols, filename, lt )
        implicit none
        integer :: file_unit, status
        integer :: nrows, ncols, i, j
        character(len=100) :: filename
        character(len=4) :: lt

        ! filename = "./models/nn-02/w-b/fc1w.txt"
        ! nrows = 8
        ! ncols = 20
        open(unit=file_unit, file=trim(filename), status='old', action='read', iostat=status)
        if (status /= 0) then
            print*, "Error opening file ", trim(filename)
            stop
        endif

        print *, "LT:"
        print *, trim(lt)
        print *, "\n"
        if (trim(lt) == "fc1w") then
            do i = 1, nrows
                read(file_unit, *) (fc1w(i, j), j = 1, ncols)
            end do
        else if (trim(lt) == "fc1b") then
            do i = 1, nrows
                read(file_unit, *) (fc1b(i, j), j = 1, ncols)
            end do
        else if (trim(lt) == "fc2w") then
            do i = 1, nrows
                read(file_unit, *) (fc2w(i, j), j = 1, ncols)
            end do
        else if (trim(lt) == "fc2b") then
            do i = 1, nrows
                read(file_unit, *) (fc2b(i, j), j = 1, ncols)
             end do
        end if

        close(file_unit)

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
        real(kind=8), intent(inout), dimension(20, 1) :: input
        integer, intent(out) :: classification

        real(kind=8), dimension(20, 20) :: fc1_w
        real(kind=8), dimension(20, 1) :: fc1_b

        real(kind=8), dimension(2, 20) :: fc2_w
        real(kind=8), dimension(2, 1) :: fc2_b

        real(kind=8), dimension(20, 1) :: temp1
        real(kind=8), dimension(2, 1) :: temp2

        real(kind=8) :: min_val, max_val

        !min_val = minval(input(1:7,1))
        !max_val = maxval(input(1:7,1))

        !input(1:7,1) = -1.0 + ( input(1:7,1) -  min_val) * 2.0 / (max_val - min_val)

        !fc1_w = reshape([-0.3494586, 1.4086386, -1.5963485, 1.1167134, -0.38775367, -0.85556597, 1.1059996, -0.671563, 0.050759546, -0.650722, 0.7240435, -1.2372845, -0.79383796, 0.5529582, 1.4097283, -1.1833408, 0.3239814, -0.31454444, 1.0651594, -0.25041237, -0.22440208, 0.9360802, -0.74616575, -0.19450942, 0.043857932, 0.8511929, 0.2386647, 1.1455798, 0.9548695, -0.20641705, -0.06948202, 0.101596326, 0.6951378, 0.42475352, 1.6859119, -0.9241148, -1.0657426, 0.19718482, 0.4095433, -1.416926, -0.012679624, 0.104560524, 0.62450576, 1.6349591, 0.9851863, -0.23168916, -0.9781, 0.8554866, -0.436441, 2.289598, 0.25891197, -0.82032084, -1.1888598, 0.049802013, 2.175047, -1.4883752, 0.32282016, -0.43852875, 0.39626214, -0.34334773, -0.23446772, -0.5600062, -0.19601557, -0.3891512, -1.2483301, 1.2013749, -1.2634141, -0.43017882, 0.30892628, 0.4304394, 1.0846772, 0.32481617, 0.96098846, 0.8961746, -0.67440623, -1.7339851, -0.38037562, -0.6060334, -0.06284777, -1.512661, -0.944193, 0.546659, 0.4180359, -0.2695023, 1.1251086, 0.50272495, 0.5559879, 0.2925843, 0.5549273, 0.6471844, -0.8365946, -0.33894426, -0.0052949768, -0.43233672, 0.22085479, -0.8741017, -0.32224286, 0.65086174, 0.53571606, 0.059319247, -1.4321054, -0.96640384, 0.072983414, -1.5505593, -1.3618819, -1.028095, -0.6055051, 0.8877591, -0.19339003, 0.2216086, 0.3123703, 0.9704603, -0.8479542, 0.35984316, -0.0042555607, 0.77339685, -0.5319276, 0.4746232, 0.47084263, -0.33044294, 0.09058956, 0.101701386, -0.56095105, -0.86824423, 0.62526786, 0.63251626, 0.38875052, -0.23737863], shape = [16,8])
        !fc1_b = reshape([-0.7604262, 0.0775104, 0.0019669298, 0.49476877, -0.0036458296, 0.16498171, 0.06969554, -0.0032051955, 0.074481755, -0.034850843, -0.4158969, -0.0030337803, 0.8863704, 0.42020145, -0.00075900427, 0.0069708955], shape = [16,1])
        !fc1_b = reshape([-0.3907799, 0.06023793, -0.36605054, -0.046650548, 0.69907707, -0.073601685, -0.50154054, 0.71009886],shape=[8,1])

        !fc2_w = reshape([-1.2119217, 1.1785002, 0.096056506, -0.09621019, -1.0274725, 0.43894228, -0.38164598, 0.78386766, 0.44864503, -1.0187023, -0.1834192, 0.19408536, -0.9547796, 0.9513067, 1.0542917, -0.72312576], shape=[2,8])

        !fc2_b = reshape([0.5561742, -0.556224], shape=[2, 1])

        temp1 = fc1b + matmul(fc1w, input)
        call l_relu(temp1, 20)
        temp2 = fc2b + matmul(fc2w, temp1)


        temp2(:, 1) = temp2(:, 1) - maxval(temp2(:, 1))
        if (exp(temp2(2, 1))/(exp(temp2(1, 1)) + exp(temp2(2, 1))) >= 0.99) then
        !if (temp2(1, 1) >= temp2(2, 1)) then
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
