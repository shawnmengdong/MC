program  MC
    USE readMatrix ! module 
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
    implicit none      ! Used for MPI stuff
    integer ArraySize,N_index,i,pivot_row,row,col,t1,t2,t3,t4,t_cal,cr
    real(kind=8), dimension(:), ALLOCATABLE :: A(:,:)     ! Create dynamic array for 16x16 matrix
    real(kind=8), dimension(:), ALLOCATABLE :: column_array(:)    ! Create dynamic array for column data communication block
    real(kind=8), dimension(:), ALLOCATABLE :: row_array(:)    ! Create dynamic array for data row communication block
    real(kind=8) :: global_log_det,pivot_value
    ! for command line arguments
    character (len=10) :: arg
    character (len=50) :: fname

    global_log_det = 0.0d+0
    t_cal = 0

    ! Get matrix size from command line
    !------------------------------------------------------------------------------------------------
    if (command_argument_count() > 0) then
            call get_command_argument(1,arg)
            READ (arg(:),'(I10)') ArraySize
    else
            print *,'Need matrix size as first parameter'
            call EXIT(1)
    endif
    !------------------------------------------------------------------------------------------------
    N_index = ArraySize
    ALLOCATE(A(ArraySize,ArraySize),column_array(ArraySize),row_array(ArraySize))  ! Allocate memory now
    ! If file name is given on command line use it
    if (command_argument_count() > 1) then
            call get_command_argument(2,fname)
            call readM(fname,ArraySize,A)
    else
    ! Otherwise fill matrix with rand values from -0.5 to 0.5
            call getRandM(ArraySize,ArraySize,A)
    endif

    !Initialize clock
    call system_clock(count_rate = cr)
    call system_clock(t1)
    !------------------------------------------------------------------------------------------------
    !This part is matrix condensation on one processor
    do i=1,N_index-2
        N_index= ArraySize-i  !  Number of rows for condensed matrix
        pivot_row = maxloc(abs(A(1:N_index+1,i)),DIM=1)  ! Find the maximum from entire column
        pivot_value = A(pivot_row,i)
        column_array(1:N_index+1) = A(1:N_index+1,i)/pivot_value !  Rest of the column is divided by pivot and stored in column array
        column_array(pivot_row) = column_array(N_index+1)   ! Switch last effective row and pivot row for column_array
        row_array(i+1:ArraySize) = A(pivot_row,i+1:ArraySize)  ! Rest of the row is stored in the row array
        global_log_det = global_log_det + log10(abs(pivot_value))  !Update global determinant
        A(pivot_row,i+1:ArraySize) = A(N_index+1,i+1:ArraySize)  !Switch last effective row and pivot row for A
        call system_clock(t2)
        do col = i+1,ArraySize
            do row = 1,N_index
                A(row,col)= A(row,col) - column_array(row)*row_array(col)
            end do
        end do
        call system_clock(t3)
        t_cal = t_cal + t3-t2
    end do
    !-------------------------------------------------------------------
    !Finally, there is only a 2 by 2 matrix left stored at A
    global_log_det = global_log_det + log10(abs(A(1,i)*A(2,i+1)-A(1,i+1)*A(2,i)))
    call system_clock(t4)
    write(*,*) 'The log10 of determinant is : ',global_log_det
    write(*,*) '1 Processors, ','Total Time used: ',real(t4-t1)/cr,' seconds'
    write(*,*) '1 Processors, ','Total Time used in calculation: ',real(t_cal)/cr,' seconds'
    DEALLOCATE(A,column_array,row_array)  ! Allocate memory now

end program MC
