program  MC
     USE MPI      ! Used to make reading the bin file easy
     USE readMatrix ! module
     USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
     implicit none      ! Used for MPI stuff
     integer myrank,numtasks,ierr,ArraySize,ArraySizep,rem,N_row,N_col,i,j,master_p,pivot_row,row,col,col_shift
     integer num_ploop,t1,t2,t3,t4,t5,t_mpi,cr
     real(kind=8), dimension(:), ALLOCATABLE :: A(:,:)     ! Create Matrix for A
     real(kind=8), dimension(:), ALLOCATABLE :: local_A(:,:)    ! Create Matrix for local_A
     real(kind=8), dimension(:), ALLOCATABLE :: local_Results(:,:)    ! Create dynamic array for local_matrix
     real(kind=8), dimension(:), target,ALLOCATABLE :: column_array(:)    ! Create dynamic array for column data communication block
     real(kind=8), dimension(:), ALLOCATABLE :: row_array(:)    ! Create dynamic array for data row communication block
     real(kind=8) :: local_log_det,global_log_det,pivot_value,sub_pivot
    ! for command line arguments
    character (len=10) :: arg
    character (len=50) :: fname

    ! Get matrix size from command line
    !------------------------------------------------------------------------------------------------
    if (command_argument_count() > 0) then
            call get_command_argument(1,arg)
            READ (arg(:),'(I10)') ArraySize
    else
            print *,'Need matrix size as first parameter'
            call EXIT(1)
    endif
     call MPI_INIT( ierr )
     call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
     local_log_det = 0.0
     global_log_det = 0.0
     t_mpi=0
     ArraySizep = ArraySize/numtasks
     N_row = ArraySize
     N_col = ArraySizep
     rem = mod(ArraySize,numtasks)
     num_ploop = numtasks-1  ! the last processor inedx is numtasks-1
     if (myrank==0) then
        N_col = N_col+rem
     endif
     !------------------------------------------------------------------------------------------------------------
     ALLOCATE(A(ArraySize,ArraySize),local_A(N_row,N_col),column_array(N_row),row_array(N_col),&
     local_Results(numtasks+rem,numtasks+rem))
     if (myrank==0) THEN
     ! If file name is given on command line use it
         if (command_argument_count() > 1) then
            call get_command_argument(2,fname)
            call readM(fname,ArraySize,A)
         else
     ! Otherwise fill matrix with rand values from -0.5 to 0.5
            call getRandM(ArraySize,ArraySize,A)
         endif
     ENDIF
     call system_clock(count_rate = cr)
     call system_clock(t1)
     !-----------------------------------------------------------------------------------------------------------
     !distribution of divisible columns
     call MPI_SCATTER(A(:,1:numtasks*ArraySizep),ArraySize*ArraySizep,MPI_DOUBLE,local_A,ArraySize*ArraySizep,&
    MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)   ! Scatter the vector to local array
     !distribution of extra columns
     if (rem>0 .And. myrank ==0) then
        local_A(:,ArraySizep+1:N_col)= A(:,numtasks*ArraySizep+1:ArraySize)  !Giving all extra to p0
     endif
     call system_clock(t5)
     !------------------------------------------------------------------------------------------------
      ! Begin local matrix condensation
      do i=1,ArraySizep-1
        do master_p = 0,num_ploop
            if (myrank==master_p) THEN
                pivot_row = maxloc(abs(local_A(1:N_row,i)),DIM=1)
                pivot_value = local_A(pivot_row,i)
                column_array(1:N_row) = local_A(1:N_row,i)/pivot_value !  Rest of the array is the normalized row
                column_array(pivot_row) = column_array(N_row)  ! Switch last effective row and pivot row for column_array
                column_array(N_row) = pivot_row  ! last element of Array is the row_index
                local_log_det = local_log_det + log10(abs(pivot_value))  !Update local determinant
                col_shift = 1
            else
                col_shift = 0
            ENDIF
            call system_clock(t2)
            call MPI_Bcast(column_array,N_row,MPI_DOUBLE,master_p,MPI_COMM_WORLD,ierr)  ! Broadcast the array to all processors
            call system_clock(t3)
            t_mpi = t_mpi+t3-t2
            !--------------------------------------------------------------------
            ! This part is matrix condensation calculation
            pivot_row = INT(column_array(N_row))
            N_row = N_row-1   ! Number of row will be one less
            row_array(i+col_shift:N_col) = local_A(pivot_row,i+col_shift:N_col)  !Get the row_array
            local_A(pivot_row,i+col_shift:N_col) = local_A(N_row+1,i+col_shift:N_col)  !Switch last effective row and pivot row for local_A
            do col = i+col_shift,N_col
                do row = 1,N_row
                     local_A(row,col)= local_A(row,col) - column_array(row)*row_array(col)
                end do
            end do
        end do
      end do
      !First copy remainder to local results on p0
      if (myrank==0 .And. rem>0) then
        local_Results(:,1:rem) = local_A(1:N_row,ArraySizep:N_col-1)
      endif
      !------------------------------------------------------------------------------------------------
      !This part is Gathering and Reduction of results
      call system_clock(t2)
      call MPI_Gather(local_A(1:N_row,N_col),N_row,MPI_DOUBLE,local_Results(:,rem+1:numtasks+rem),N_row,MPI_DOUBLE,0,&
      MPI_COMM_WORLD,ierr)
      call MPI_Reduce(local_log_det,global_log_det,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)  ! Reduce determinants to master nodes
      call system_clock(t3)
      t_mpi = t_mpi+t3-t2
      !------------------------------------------------------------------------------------------------
      !This part is matrix condensation on one processor
       if (myrank==0) THEN
         do i=1,numtasks+rem
             pivot_value = local_Results(i,i)
             global_log_det = global_log_det + log10(abs(pivot_value))  !Update global determinant
             column_array(i+1:numtasks+rem) = local_Results(i+1:numtasks+rem,i)/pivot_value   !Normalize the column
             do j = i+1,numtasks+rem
                 sub_pivot = local_Results(i,j)
                 local_Results(i+1:numtasks+rem,j)= local_Results(i+1:numtasks+rem,j)-sub_pivot*column_array(i+1:numtasks+rem)
             end do
         end do
         call system_clock(t4)
         write(*,*) 'The log10 of determinant is : ',global_log_det
         write(*,*) numtasks,'Processors, ','Total Time used: ',real(t4-t1)/cr,' seconds'
         write(*,*) numtasks,'Processors, ','Total Communication Time used: ',real(t_mpi)/cr,' seconds'
         write(*,*) numtasks,'Processors, ','Data Distribution Time used: ',real(t5-t1)/cr,' seconds'
      endif
      DEALLOCATE(A,local_A,local_Results,column_array,row_array)  ! Allocate memory now
      call MPI_FINALIZE(ierr)
end program MC
