program  GE_col
    USE readMatrix ! module
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
    implicit none      ! Used for MPI stuff
    integer ArraySize,i,j,t1,t4,cr
    real(kind=8), dimension(:), ALLOCATABLE :: A(:,:)     ! Create dynamic array for 16x16 matrix
    real(kind=8), dimension(:), ALLOCATABLE :: column_array(:)
    real(kind=8) :: global_log_det,pivot,sub_pivot
    ! for command line arguments
    character (len=10) :: arg
    character (len=50) :: fname

    global_log_det = 0.0d+0

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
    ALLOCATE(A(ArraySize,ArraySize),column_array(ArraySize))  ! Allocate memory now
    ! If file name is given on command line use it
    if (command_argument_count() > 1) then
            call get_command_argument(2,fname)
            call readM(fname,ArraySize,A)
    else
    ! Otherwise fill matrix with rand values from -0.5 to 0.5
            call getRandM(ArraySize,ArraySize,A)
    endif
    !------------------------------------------------------------------------------------------------
    !Initialize clock
    call system_clock(count_rate = cr)
    call system_clock(t1)
    !------------------------------------------------------------------------------------------------
    !This part is matrix condensation on one processor
    do i=1,ArraySize
        pivot = A(i,i)
        global_log_det = global_log_det + log10(abs(pivot))  !Update global determinant
        column_array(i+1:ArraySize) = A(i+1:ArraySize,i)/pivot   !Normalize the column
        do j = i+1,ArraySize
            sub_pivot = A(i,j)
            A(i+1:ArraySize,j)= A(i+1:ArraySize,j) - sub_pivot*column_array(i+1:ArraySize)
        end do
    end do
    !-------------------------------------------------------------------
    call system_clock(t4)
    write(*,*) 'The log10 of determinant is : ',global_log_det
    write(*,*) '1 Processors, ','Total Time used: ',real(t4-t1)/cr,' seconds'
    DEALLOCATE(A,column_array)  ! Allocate memory now

end program GE_col
