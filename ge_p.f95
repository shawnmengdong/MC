program  GE_p
     USE readMatrix ! module
     USE MPI      ! Used to make reading the bin file easy
     USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE
     implicit none
     integer myrank,numtasks,ierr,ArraySize,ArraySizep,rem,N_row,N_col,i,master_p,pivot_column,col,col_shift
     integer num_outerloop,num_ploop,t1,t2,t3,t4,t5,t_mpi,cr,global_step
     real(kind=8), dimension(:), ALLOCATABLE :: A(:,:)     ! Create Matrix for A
     real(kind=8), dimension(:), ALLOCATABLE :: local_A(:,:)    ! Create Matrix for local_A
     real(kind=8), dimension(:), target,ALLOCATABLE :: column_array(:)    ! Create dynamic array for column data communication block
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
     global_step = 0
     ArraySizep = ArraySize/numtasks
     rem = mod(ArraySize,numtasks)
     num_ploop = numtasks-1  ! the last processor inedx is numtasks-1
     num_outerloop=ArraySizep
     !------------------------------------------------------------------------------------------------------------------
     if (myrank==0) then
        ArraySizep = ArraySizep+rem
     endif
     N_row = ArraySize
     N_col = ArraySizep

     ALLOCATE(A(ArraySize,ArraySize),local_A(N_row,N_col),column_array(N_row))  ! Allocate memory now

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
     call system_clock(t1) !Begin clock
    !--------------------------------------------------------------------------------------------------
    !Cyclic Distribution
     do i = 1,num_outerloop
         call MPI_SCATTER(A(:,numtasks*(i-1)+1:i*numtasks),ArraySize,MPI_DOUBLE,local_A(:,i),ArraySize, &
         MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
     end do
     if (rem>0 .AND. myrank ==0)then
         local_A(:,num_outerloop+1:ArraySizep)= A(:,numtasks*num_outerloop+1:ArraySize)  !Giving all extra to p0
     endif

     call system_clock(t5)
     !------------------------------------------------------------------------------------------------
     ! Begin local Gaussian Elimination
      do i=1,num_outerloop
        do master_p = 0,num_ploop
            global_step = global_step +1
            if (myrank==master_p) THEN
                pivot_column = i+maxloc(abs(local_A(global_step,i:N_col)),DIM=1)-1
                pivot_value = local_A(global_step,pivot_column)
                column_array(global_step+1:ArraySize) = local_A(global_step+1:ArraySize,pivot_column)/pivot_value !  Rest of the array is the normalized row
                local_A(global_step:ArraySize,pivot_column)=local_A(global_step:ArraySize,i)  !Switch columns
                local_log_det = local_log_det + log10(abs(pivot_value))  !Update local determinant
                col_shift = 1
            else
                col_shift = 0
            ENDIF
            N_row = N_row-1 ! Number of row will be one less


            call system_clock(t2)
            call MPI_Bcast(column_array(global_step+1:ArraySize),N_row,MPI_DOUBLE,master_p,MPI_COMM_WORLD,ierr)  ! Broadcast the array to all processors
            call system_clock(t3)
            t_mpi = t_mpi+t3-t2
            !--------------------------------------------------------------------
            ! This part is local gaussian elimination
            do col = i+col_shift,N_col
                sub_pivot = local_A(global_step,col)
                local_A(global_step+1:ArraySize,col) = local_A(global_step+1:ArraySize,col) - &
                sub_pivot * column_array(global_step+1:ArraySize)
            end do
        end do
      end do
      call system_clock(t2)
      call MPI_Reduce(local_log_det,global_log_det,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)  ! Reduce determinants to master nodes
      call system_clock(t3)
      t_mpi = t_mpi+t3-t2
      !------------------------------------------------------------------------------------------------
      !This part is Gaussian Elimination on one processor
      if (myrank==0) THEN
         do i = num_outerloop+1,ArraySizep-1
            global_step = global_step +1
            !partial pivoting:
            pivot_column = i+maxloc(abs(local_A(global_step,i:N_col)),DIM=1)-1
            pivot_value = local_A(global_step,pivot_column)
            column_array(global_step+1:ArraySize) = local_A(global_step+1:ArraySize,pivot_column)/pivot_value !  Rest of the array is the normalized row
            local_A(global_step:ArraySize,pivot_column)=local_A(global_step:ArraySize,i)  !Switch columns
            global_log_det = global_log_det + log10(abs(pivot_value))  !Update local determinant

            !elimination
            do col = i,N_col
                sub_pivot = local_A(global_step,col)
                local_A(global_step+1:ArraySize,col) = local_A(global_step+1:ArraySize,col) - &
                sub_pivot * column_array(global_step+1:ArraySize)
            end do
         end do
         if (rem>0) then
            global_log_det = global_log_det + log10(abs(local_A(ArraySize,ArraySizep)))
         end if

         call system_clock(t4)
         write(*,*) 'The log10 of determinant is : ',global_log_det
         write(*,*) numtasks,'Processors, ','Total Time used: ',real(t4-t1)/cr,' seconds'
         write(*,*) numtasks,'Processors, ','Total Communication Time used: ',real(t_mpi)/cr,' seconds'
         write(*,*) numtasks,'Processors, ','Matrix Distribution time: ',real(t5-t1)/cr,' seconds'
      ENDIF
      DEALLOCATE(A,local_A,column_array)  ! Allocate memory now
      call MPI_FINALIZE(ierr)
end program GE_p
