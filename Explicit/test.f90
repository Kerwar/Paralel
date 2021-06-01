program test
  
  use type_m   , only: DP 
  use param_m  , only: param_t 
  use paralel_m, only: parll_t
  
  implicit none
  
  include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'
  
  type(parll_t) :: PL
  type(param_t) :: param
  
  integer :: ierr, i, j
  
  integer, allocatable :: i_str_p(:), i_end_p(:), &
  j_str_p(:), j_end_p(:)
  ! Initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, PL%nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, PL%myproc, ierr)
  
  call param%read_param()
  
  call PL%read_parll()
  
  ! Domain decomposition in processors
  PL%nrows = int(sqrt(float(param%ny)/float(param%nx)* &
  float(PL%nprocs)))
  if (PL%nrows == 0) PL%nrows = 1
  
  PL%ncols = int(float(PL%nprocs)/float(PL%nrows))
  
  do while (PL%nrows * PL%ncols /= PL%nprocs)
    PL%nrows = PL%nrows - 1
    PL%ncols = int(float(PL%nprocs)/float(PL%nrows))
  end do
  
  ! Domain decompositions
  
  PL%nxadd = mod(param%nx, PL%ncols)
  PL%nyadd = mod(param%ny, PL%nrows)
  
  PL%myrow = min(int(float(PL%myproc)/float(PL%ncols)), PL%ncols-1)
  PL%mycol = PL%myproc - PL%myrow * PL%ncols 
  
  ! Set size of every mesh for each processor
  if (PL%mycol < PL%nxadd) then
    PL%istr = PL%mycol * &
    (int(float(param%nx)/float(PL%ncols)) + 1) + 1
    PL%iend = (PL%mycol + 1) * &
    (int(float(param%nx)/float(PL%ncols)) + 1) 
  else
    PL%istr = PL%nxadd + PL%mycol * &
    int(float(param%nx)/float(PL%ncols))  + 1
    PL%iend = PL%nxadd + (PL%mycol + 1) * &
    int(float(param%nx)/float(PL%ncols)) + 1
  end if
  
  if (PL%myrow < PL%nyadd) then
    PL%jstr = PL%myrow * &
    (int(float(param%ny)/float(PL%nrows)) + 1) + 1
    PL%jend = (PL%myrow + 1) * &
    (int(float(param%ny)/float(PL%nrows)) + 1) 
  else
    PL%jstr = PL%nyadd + PL%myrow * &
    int(float(param%ny)/float(PL%nrows)) + 1
    PL%jend = PL%nyadd + (PL%myrow + 1) * &
    int(float(param%ny)/float(PL%nrows)) + 1 
  end if
  
  PL%nx = PL%iend - PL%istr + 1
  PL%ny = PL%jend - PL%jstr + 1
  
  print *, "My Proc:", PL%myproc, "My Row:", PL%myrow, "My Col:",PL%mycol, &
  "Istr:", PL%istr, "Iend:", PL%iend, "Jstr:",PL%jstr, "Jend:", PL%jend
  
  if (PL%myproc == 0) then
    allocate(i_str_p(0:PL%nprocs-1), i_end_p(0:PL%nprocs-1), & 
    j_str_p(0:PL%nprocs-1), j_end_p(0:PL%nprocs-1))
    
    do j = 0, PL%nrows-1
      do i = 0, PL%nxadd-1
        i_str_p(i+j*PL%ncols) =  i    * (param%nx/PL%ncols + 1) + 1 
        i_end_p(i+j*PL%ncols) = (i+1) * (param%nx/PL%ncols + 1) 
      enddo
      do i = PL%nxadd, PL%ncols-1
        i_str_p(i+j*PL%ncols) = PL%nxadd + &
         i      * (param%nx/PL%ncols) + 1
        i_end_p(i+j*PL%ncols) = PL%nxadd + &
        (i + 1) * (param%nx/PL%ncols) 
      enddo
    enddo
    
    do i = 0, PL%ncols -1
      do j = 0, PL%nyadd-1
        j_str_p(j*PL%ncols+i) =  j    * (param%ny/PL%nrows + 1) + 1
        j_end_p(j*PL%ncols+i) = (j+1) * (param%ny/PL%nrows + 1) 
      enddo
      do j = PL%nyadd, PL%nrows-1
        j_str_p(j*PL%ncols+i) = PL%nyadd + &
        j     * (param%ny/PL%nrows) + 1
        j_end_p(j*PL%ncols+i) = PL%nyadd + &
        (j+1) * (param%ny/PL%nrows)
      enddo
    enddo

    print *, j_str_p(:)
    print *, j_end_p(:)
    
    deallocate(i_str_p, i_end_p, j_str_p, j_end_p)
  end if
  
  call MPI_FINALIZE(ierr)
end program test