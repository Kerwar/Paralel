program main
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SET UP OF THE PROBLEM, need to take in account there is a simmetry.
  ! We only consider half of each channel and then we impose cubic simmetry
  !
  ! ----------------------------------------------------------------------
  ! 
  !  ---------------->
  !
  ! ----------------------------------------------------------------------
  !                         |                |
  ! ----------------------------------------------------------------------
  !
  !                                                      <----------------
  !
  ! ----------------------------------------------------------------------
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Processors numeration  
  !  _____________________________________________________________________
  ! |         6         |              7           |           8          |
  ! |___________________|__________________________|______________________|
  ! |         3         |              4           |           5          |
  ! |___________________|__________________________|______________________|
  ! |         0         |              1           |           2          |
  ! |___________________|__________________________|______________________|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use type_m   , only: DP 
  use param_m  , only: param_t 
  use paralel_m, only: parll_t
  use coefun_m , only: ADI, write_tstep
  use point_m  , only: point_t, coefx_1, coefx_2, coefx_3, coefx_4, coefx_5, coefx_6, &
  coefx_7, coefx_8, coefx_9, coefy_1, coefy_2, coefy_3, coefy_4, coefy_5, init_gamma
  
  implicit none
  
  include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'
  
  type(param_t) :: param
  type(parll_t) :: PL
  
  real(DP) :: t_cpu_start, t_cpu_end
  
  type(point_t), allocatable :: field(:,:), field0(:,:), solution(:,:,:)
  ! Message Send/Recive Left/Riht/Bot/Top
  real(DP), allocatable :: msl(:), mrl(:), msr(:), mrr(:), &
  msb(:), mrb(:), mst(:), mrt(:), soltransf(:), everyerror(:)
  
  real(DP) :: tmperror, error
  ! Error for MPI
  integer :: ierr!, ROW_COMM, COL_COMM
  
  integer :: status(MPI_STATUS_SIZE)
  
  ! Integer needed for indeces
  integer :: i, j, k, l, t, meshsize, rowsize, colsize, proc, index
  ! Indeces for data treatment in proc0
  integer, allocatable :: i_str_p(:), i_end_p(:), &
  j_str_p(:), j_end_p(:)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI SET UP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Initialize MPI
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, PL%nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, PL%myproc, ierr)
  
  if (PL%myproc == 0) call cpu_time(t_cpu_start)
  
  call param%read_param()
  
  call PL%read_parll()
  
  if (PL%myproc == 0) allocate(solution(param%nx, param%ny, param%nt), &
  everyerror(PL%nprocs))
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Domain decomposition in processors
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PL%nrows = int(sqrt(float(param%ny)/float(param%nx)* &
  float(PL%nprocs)))
  if (PL%nrows == 0) PL%nrows = 1
  
  PL%ncols = int(float(PL%nprocs)/float(PL%nrows))
  
  do while (PL%nrows * PL%ncols /= PL%nprocs)
    PL%nrows = PL%nrows - 1
    PL%ncols = int(float(PL%nprocs)/float(PL%nrows))
  end do

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Domain decomposition and setting meshing sizes
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  !NEED TO CHECK IF THIS DECOMPOSITION WORKS FOR MORE THAN 2 PROC
  PL%nxadd = mod(param%nx, PL%ncols)
  PL%nyadd = mod(param%ny, PL%nrows)
  
  PL%myrow = min(int(float(PL%myproc)/float(PL%ncols)), PL%ncols-1)
  PL%mycol = PL%myproc - PL%myrow * PL%ncols 
  
  if (PL%mycol < PL%nxadd) then
    PL%istr = PL%mycol * &
    (int(float(param%nx)/float(PL%ncols)) + 1) + 1
    PL%iend = (PL%mycol + 1) * &
    (int(float(param%nx)/float(PL%ncols)) + 1) 
  else
    PL%istr = PL%nxadd + PL%mycol * &
    int(float(param%nx)/float(PL%ncols))  + 1
    PL%iend = PL%nxadd + (PL%mycol + 1) * &
    int(float(param%nx)/float(PL%ncols)) 
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
    int(float(param%ny)/float(PL%nrows)) 
  end if
  
  PL%nx = PL%iend - PL%istr + 1
  PL%ny = PL%jend - PL%jstr + 1
  
  allocate(soltransf(3 * PL%nx * PL%ny))
  
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
  end if
  !print *, PL%myproc, PL%myrow, PL%nrows, PL%mycol, &
  !PL%ncols, PL%istr, PL%iend, PL%jstr, PL%jend
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Seting who are the neighbours
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! call MPI_COMM_SPLIT(&
  ! MPI_COMM_WORLD, PL%myrow, PL%mycol, ROW_COMM, ierr)
  ! call MPI_COMM_RANK (ROW_COMM, PL%myid_row , ierr)
  ! call MPI_COMM_SIZE (ROW_COMM, PL%nodes_row, ierr)
  
  ! call MPI_COMM_SPLIT(&
  ! MPI_COMM_WORLD, PL%mycol, PL%myrow, COL_COMM, ierr)
  ! call MPI_COMM_RANK (COL_COMM, PL%myid_col , ierr)
  ! call MPI_COMM_SIZE (COL_COMM, PL%nodes_col, ierr)
  
  !Set IDS of neighbours
  PL%myleft = PL%myproc - 1
  if (mod(PL%myproc , PL%ncols) == 0) PL%myleft  = MPI_PROC_NULL
  
  PL%myright = PL%myproc + 1
  if (mod(PL%myright, PL%ncols) == 0) PL%myright = MPI_PROC_NULL
  
  ! We impose simmetry
  if (PL%nrows /= 1) then
    PL%mybot = PL%myproc - PL%ncols
    if (PL%mybot < 0) PL%mybot = PL%nprocs - PL%myproc - 1
    
    PL%mytop = PL%myproc + PL%ncols
    if (PL%mytop >= PL%nprocs) PL%mytop = PL%ncols - mod(PL%myproc, PL%ncols) - 1
  else 
    PL%mybot = MPI_PROC_NULL
    
    PL%mytop = MPI_PROC_NULL
  end if
  
  if (PL%myproc == 0) call cpu_time(t_cpu_end)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialization of the field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  allocate(field(0:PL%nx+1, 0:PL%ny+1))
  allocate(field0(0:PL%nx+1, 0:PL%ny+1))
  
  !print *, "Memory allocated in the ", PL%myproc, " processor!"

    field(0,0)%aT(:) = 0.0_DP
    field(0,0)%aF(:) = 0.0_DP
    field(0,0)%aZ(:) = 0.0_DP
    
    field(0,0)%aT(5) = 1.0_DP
    field(0,0)%aF(5) = 1.0_DP
    field(0,0)%aZ(5) = 1.0_DP

    field(0,0)%sT = 0.0_DP
    field(0,0)%sF = 0.0_DP
    field(0,0)%sZ = 0.0_DP

    field(PL%nx + 1,0)%aT(:) = field(0,0)%aT(:)
    field(PL%nx + 1,0)%aF(:) = field(0,0)%aF(:)
    field(PL%nx + 1,0)%aZ(:) = field(0,0)%aZ(:)

    field(PL%nx + 1,0)%sT = 0.0_DP
    field(PL%nx + 1,0)%sF = 0.0_DP
    field(PL%nx + 1,0)%sZ = 0.0_DP
  do i = 1, PL%nx

    field(i,0)%aT(:) = field(0,0)%aT(:)
    field(i,0)%aF(:) = field(0,0)%aF(:)
    field(i,0)%aZ(:) = field(0,0)%aZ(:)

    do j = 1, PL%ny
      field(0,j)%aT(:) = field(0,0)%aT(:)
      field(0,j)%aF(:) = field(0,0)%aF(:)
      field(0,j)%aZ(:) = field(0,0)%aZ(:)

      field(0,j)%sT = 0.0_DP
      field(0,j)%sF = 1.0_DP
      field(0,j)%sZ = 0.0_DP
      
      call field(i,j)%set_coef(param, i + PL%istr - 1, &
      j + PL%jstr - 1)
      
      call init_gamma(field(i,j), param)
      !if (PL%myproc == 0) print *, i,j, field(i,j)%T, field(i,j)%F, field(i,j)%Z
      
      field(PL%nx + 1,j)%aT(:) = field(0,0)%aT(:)
      field(PL%nx + 1,j)%aF(:) = field(0,0)%aF(:)
      field(PL%nx + 1,j)%aZ(:) = field(0,0)%aZ(:)
      
      field(PL%nx + 1,j)%sT = 0.0_DP
      field(PL%nx + 1,j)%sF = 1.0_DP
      field(PL%nx + 1,j)%sZ = 0.0_DP
    end do
    
    field(i,PL%ny + 1)%aT(:) = field(0,0)%aT(:)
    field(i,PL%ny + 1)%aF(:) = field(0,0)%aF(:)
    field(i,PL%ny + 1)%aZ(:) = field(0,0)%aZ(:)
    
  end do
  
    field(0        ,PL%ny + 1)%aT(:) = field(0,0)%aT(:)
    field(PL%nx + 1,PL%ny + 1)%aT(:) = field(0,0)%aT(:)
    field(0        ,PL%ny + 1)%aF(:) = field(0,0)%aF(:)
    field(PL%nx + 1,PL%ny + 1)%aF(:) = field(0,0)%aF(:)
    field(0        ,PL%ny + 1)%aZ(:) = field(0,0)%aZ(:)
    field(PL%nx + 1,PL%ny + 1)%aZ(:) = field(0,0)%aZ(:)
  
  field0 = field

  allocate(msl(3 * PL%ny), mrl(3 * PL%ny))
  allocate(msr(3 * PL%ny), mrr(3 * PL%ny))
  allocate(msb(3 * PL%nx), mrb(3 * PL%nx))
  allocate(mst(3 * PL%nx), mrt(3 * PL%nx))
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Sending info to neightbours
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 0,1
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info left
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_col, 2) == i) then
      do j = 1, PL%ny
        msl(j)             = field(1,j)%T
        msl(j + PL%ny)     = field(1,j)%F
        msl(j + 2 * PL%ny) = field(1,j)%Z
      end do
      
      call MPI_SEND(msl, size(msl), MPI_DOUBLE_PRECISION, &
      PL%myleft, 0, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrr, size(mrr), MPI_DOUBLE_PRECISION, &
      PL%myright, 0, MPI_COMM_WORLD, status,  ierr)
      
      do j = 1, PL%ny
        field(PL%nx + 1,j)%T  = mrr(j)
        field(PL%nx + 1,j)%sT = field(PL%nx + 1,j)%T

        field(PL%nx + 1,j)%F  = mrr(j + PL%ny)
        field(PL%nx + 1,j)%sF = field(PL%nx + 1,j)%F

        field(PL%nx + 1,j)%Z  = mrr(j + 2 * PL%ny)
        field(PL%nx + 1,j)%sZ = field(PL%nx + 1,j)%Z

      end do
    end if
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info right
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_col, 2) == i) then
      do j = 1, PL%ny
        msr(j)             = field(PL%nx,j)%T
        msr(j + PL%ny)     = field(PL%nx,j)%F
        msr(j + 2 * PL%ny) = field(PL%nx,j)%Z
      end do
        
      call MPI_SEND(msr, size(msr), MPI_DOUBLE_PRECISION, &
      PL%myright, 1, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrl, size(mrl), MPI_DOUBLE_PRECISION, &
      PL%myleft, 1, MPI_COMM_WORLD, status, ierr)
      
      do j = 1, PL%ny
        field(0,j)%T = mrl(j)
        field(0,j)%sT = field(0,j)%T

        field(0,j)%F = mrl(j + PL%ny)
        field(0,j)%sF = field(0,j)%F

        field(0,j)%Z = mrl(j + 2* PL%ny)
        field(0,j)%sZ = field(0,j)%Z
      end do
    end if
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info bot
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_row, 2) == i) then
      do j = 1, PL%nx
        msb(j)             = field(j,1)%T
        msb(j + PL%nx)     = field(j,1)%F
        msb(j + 2 * PL%nx) = field(j,1)%Z
      end do
        
      call MPI_SEND(msb, size(msb), MPI_DOUBLE_PRECISION, &
      PL%mybot, 2, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrt, size(mrt), MPI_DOUBLE_PRECISION, &
      PL%mytop, 2, MPI_COMM_WORLD, status, ierr)
        
      do j = 1, PL%nx
        field(j,PL%ny + 1)%T  = mrt(j)
        if(PL%myrow == 0) field(j,PL%ny + 1)%T = mrt(PL%nx - j + 1)
        field(j,PL%ny + 1)%sT = field(j,PL%ny + 1)%T

        field(j,PL%ny + 1)%F  = mrt(j + PL%nx)
        if(PL%myrow == 0) field(j,PL%ny + 1)%F = mrt(2 * PL%nx - j + 1)
        field(j,PL%ny + 1)%sF = field(j,PL%ny + 1)%F

        field(j,PL%ny + 1)%Z  = mrt(j + 2 * PL%nx)
        if(PL%myrow == 0) field(j,PL%ny + 1)%Z = mrt(3 * PL%nx - j + 1)
        field(j,PL%ny + 1)%sZ = field(j,PL%ny + 1)%Z

      end do
    end if
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info top
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_row, 2) == i) then
      do j = 1, PL%nx
        mst(j)             = field(j,PL%ny)%T
        mst(j + PL%nx)     = field(j,PL%ny)%F
        mst(j + 2 * PL%nx) = field(j,PL%ny)%Z
      end do
        
      call MPI_SEND(mst, size(mst), MPI_DOUBLE_PRECISION, &
      PL%mytop, 3, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrb, size(mrb), MPI_DOUBLE_PRECISION, &
      PL%mybot, 3, MPI_COMM_WORLD, status, ierr)
      do j = 1, PL%nx
        field(j,0)%T  = mrb(j)
        if(PL%myrow == PL%ncols - 1) field(j, 0)%T = mrb(PL%nx - j + 1)
        field(j,0)%sT = field(j,0)%T

        field(j,0)%F  = mrb(j + PL%nx)
        if(PL%myrow == PL%ncols - 1) field(j, 0)%F = mrb(2 * PL%nx - j + 1)
        field(j,0)%sF = field(j,0)%F

        field(j,0)%Z  = mrb(j + 2 * PL%nx)
        if(PL%myrow == PL%ncols - 1) field(j, 0)%Z = mrb(3 * PL%nx - j + 1)
        field(j,0)%sZ = field(j,0)%Z

      end do
    end if
      
  end do
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Sending info to proc 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (PL%myproc /= 0) then
    !print *, "Processor ", PL%myproc, " is up to duty"
    do j = 1, PL%ny
      do i = 1, PL%nx
        index = (j - 1) * PL%nx + i
        soltransf(                    index) = field(i,j)%T
        soltransf(    PL%nx * PL%ny + index) = field(i,j)%F
        soltransf(2 * PL%nx * PL%ny + index) = field(i,j)%Z
      !  if (field(i,j)%T == 0.0_DP .and. (j > param%my2 .or. j < param%my1)) print*, i, j
      end do
    end do
    
    print *, "Actual size of ", PL%myproc, ":", size(soltransf)
    !print *, "Sending info to main proc"
    !print *, soltransf
    call MPI_SEND(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
    0, 1000 + PL%myproc, MPI_COMM_WORLD, ierr)
    !print *, "Proc" , PL%myproc, "done sending info to main prco"
  else 
    !print *, "Processor 0 is up to duty"
    t = 1
    do i = 0, PL%nrows-1
      do j = 0, PL%ncols-1
        
        proc = i * PL%nrows + j
        !print *, "We are in proc", proc, "if u ask proc 0", PL%nrows, PL%ncols
        
        rowsize = i_end_p(proc) - i_str_p(proc) + 1
        colsize = j_end_p(proc) - j_str_p(proc) + 1
        meshsize = rowsize * colsize
        
        if (i == 0 .and. j == 0) then 
          do k = 1,rowsize
            do l = 1,colsize
              
              solution(k ,l, t)%T = field(k,l)%T
              solution(k ,l, t)%F = field(k,l)%F
              solution(k ,l, t)%Z = field(k,l)%Z
            end do
          end do
          
        else
          
          deallocate(soltransf)
          allocate(soltransf( 3 * meshsize))
          print *, "Main proc idea of size of proc ", proc, ":", size(soltransf) 
          call MPI_RECV(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
          proc, 1000 + proc, MPI_COMM_WORLD, status, ierr)            
          !print *, "We are chilling in proc 0"
          do l = 1,colsize
            do k = 1,rowsize
              index = (l - 1) * rowsize + k
              
              solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%T = &
              soltransf(               index)
              solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%F = &
              soltransf(    meshsize + index)
              solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%Z = &
              soltransf(2 * meshsize + index)
            end do
          end do
        end if
        
        
      end do
    end do
  end if

  if (PL%myproc == 0) call write_tstep(solution(:,:,t), param, 0)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Main Loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do t = 1, param%nt
    
    !print *, "Processor:", PL%myproc, param%nt, PL%nx, PL%ny, param%my1, param%my2
    do i = 1, PL%nx
      do j = 1, PL%ny
        if (.not. associated(field(i,j)%setXc)) print *, i, j
        call field(i,j)%setXc(param)
        call field(i,j)%setYc(param)

        field(i,j)%sT = field(i,j)%sT + field0(i,j)%T * param%txy
        field(i,j)%sF = field(i,j)%sF + field0(i,j)%F * param%txy
        field(i,j)%sZ = field(i,j)%sZ + field0(i,j)%Z * param%txy
      end do
    end do
    
    !print *, "Coeficien functions done in the ", PL%myproc, " processor!"

    !  if (PL%myproc == 0 ) then
    !      open(102, file = "testing.txt", position = "append") 
    !     ! if (associated(field(i,j)%setXc, coefx_1)) write(102,*) "1x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_2)) write(102,*) "2x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_3)) write(102,*) "3x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_4) )write(102,*) "4x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_5)) write(102,*) "5x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_6)) write(102,*) "6x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_7)) write(102,*) "7x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_8)) write(102,*) "8x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setXc, coefx_9)) write(102,*) "9x", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setYc, coefy_1)) write(102,*) "1y", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setYc, coefy_2)) write(102,*) "2y", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setYc, coefy_3)) write(102,*) "3y", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setYc, coefy_4)) write(102,*) "4y", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     ! if (associated(field(i,j)%setYc, coefy_5)) write(102,*) "5y", i + PL%istr - 1, &
    !     ! j + PL%jstr - 1
    !     do i = 0,PL%nx + 1
    !     do j = 0,PL%ny + 1
    !     write(102,"(2I4, 6ES20.8)") i, j, field(i,j)%aZ(:), field(i,j)%sZ
    !     end do 
    !     end do
    !     close(102)
    !    end if
    ! stop
    call ADI(field(0:PL%nx+1,0:PL%ny+1)%aT(1), field(0:PL%nx+1,0:PL%ny+1)%aT(2), field(0:PL%nx+1,0:PL%ny+1)%aT(3), &
    field(0:PL%nx+1,0:PL%ny+1)%aT(4), field(0:PL%nx+1,0:PL%ny+1)%aT(5), field(:,:)%sT, field(:,:)%T, &
    tmperror, 4, param%valpha)
    error = max(tmperror, 0.0_DP)
    if(mod(t,param%show) == 0) print *, "T System solved in the ", PL%myproc, " processor!, ", error
    
    
    call ADI(field(0:PL%nx+1,0:PL%ny+1)%aF(1), field(0:PL%nx+1,0:PL%ny+1)%aF(2), field(0:PL%nx+1,0:PL%ny+1)%aF(3), &
    field(0:PL%nx+1,0:PL%ny+1)%aF(4), field(0:PL%nx+1,0:PL%ny+1)%aF(5), field(:,:)%sF, field(:,:)%F, &
    tmperror, 4, param%valpha)
    error = max(tmperror, error)
    if(mod(t,param%show) == 0)  print *, "F System solved in the ", PL%myproc, " processor!, ", error

    call ADI(field(0:PL%nx+1,0:PL%ny+1)%aZ(1), field(0:PL%nx+1,0:PL%ny+1)%aZ(2), field(0:PL%nx+1,0:PL%ny+1)%aZ(3), &
    field(0:PL%nx+1,0:PL%ny+1)%aZ(4), field(0:PL%nx+1,0:PL%ny+1)%aZ(5), field(:,:)%sZ, field(:,:)%Z, &
    tmperror, 4, param%valpha)
    error = max(tmperror, error)
    
    if(mod(t,param%show) == 0) print *, "Z System solved in the ", PL%myproc, " processor!, ", error
    
  do i = 0,1
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info left
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_col, 2) == i) then
      do j = 1, PL%ny
        msl(j)             = field(1,j)%T
        msl(j + PL%ny)     = field(1,j)%F
        msl(j + 2 * PL%ny) = field(1,j)%Z
      end do
      
      call MPI_SEND(msl, size(msl), MPI_DOUBLE_PRECISION, &
      PL%myleft, 0, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrr, size(mrr), MPI_DOUBLE_PRECISION, &
      PL%myright, 0, MPI_COMM_WORLD, status,  ierr)
      
      do j = 1, PL%ny
        field(PL%nx + 1,j)%T  = mrr(j)
        field(PL%nx + 1,j)%sT = field(PL%nx + 1,j)%T

        field(PL%nx + 1,j)%F  = mrr(j + PL%ny)
        field(PL%nx + 1,j)%sF = field(PL%nx + 1,j)%F

        field(PL%nx + 1,j)%Z  = mrr(j + 2 * PL%ny)
        field(PL%nx + 1,j)%sZ = field(PL%nx + 1,j)%Z

      end do
    end if
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info right
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_col, 2) == i) then
      do j = 1, PL%ny
        msr(j)             = field(PL%nx,j)%T
        msr(j + PL%ny)     = field(PL%nx,j)%F
        msr(j + 2 * PL%ny) = field(PL%nx,j)%Z
      end do
        
      call MPI_SEND(msr, size(msr), MPI_DOUBLE_PRECISION, &
      PL%myright, 1, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrl, size(mrl), MPI_DOUBLE_PRECISION, &
      PL%myleft, 1, MPI_COMM_WORLD, status, ierr)
      
      do j = 1, PL%ny
        field(0,j)%T = mrl(j)
        field(0,j)%sT = field(0,j)%T

        field(0,j)%F = mrl(j + PL%ny)
        field(0,j)%sF = field(0,j)%F

        field(0,j)%Z = mrl(j + 2* PL%ny)
        field(0,j)%sZ = field(0,j)%Z
      end do
    end if
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info bot
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_row, 2) == i) then
      do j = 1, PL%nx
        msb(j)             = field(j,1)%T
        msb(j + PL%nx)     = field(j,1)%F
        msb(j + 2 * PL%nx) = field(j,1)%Z
      end do
        
      call MPI_SEND(msb, size(msb), MPI_DOUBLE_PRECISION, &
      PL%mybot, 2, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrt, size(mrt), MPI_DOUBLE_PRECISION, &
      PL%mytop, 2, MPI_COMM_WORLD, status, ierr)
        
      do j = 1, PL%nx
        field(j,PL%ny + 1)%T  = mrt(j)
        if(PL%myrow == 0) field(j,PL%ny + 1)%T = mrt(PL%nx - j + 1)
        field(j,PL%ny + 1)%sT = field(j,PL%ny + 1)%T

        field(j,PL%ny + 1)%F  = mrt(j + PL%nx)
        if(PL%myrow == 0) field(j,PL%ny + 1)%F = mrt(2 * PL%nx - j + 1)
        field(j,PL%ny + 1)%sF = field(j,PL%ny + 1)%F

        field(j,PL%ny + 1)%Z  = mrt(j + 2 * PL%nx)
        if(PL%myrow == 0) field(j,PL%ny + 1)%Z = mrt(3 * PL%nx - j + 1)
        field(j,PL%ny + 1)%sZ = field(j,PL%ny + 1)%Z

      end do
    end if
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Sending info top
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mod(PL%myid_row, 2) == i) then
      do j = 1, PL%nx
        mst(j)             = field(j,PL%ny)%T
        mst(j + PL%nx)     = field(j,PL%ny)%F
        mst(j + 2 * PL%nx) = field(j,PL%ny)%Z
      end do
        
      call MPI_SEND(mst, size(mst), MPI_DOUBLE_PRECISION, &
      PL%mytop, 3, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mrb, size(mrb), MPI_DOUBLE_PRECISION, &
      PL%mybot, 3, MPI_COMM_WORLD, status, ierr)
      do j = 1, PL%nx
        field(j,0)%T  = mrb(j)
        if(PL%myrow == PL%ncols - 1) field(j, 0)%T = mrb(PL%nx - j + 1)
        field(j,0)%sT = field(j,0)%T

        field(j,0)%F  = mrb(j + PL%nx)
        if(PL%myrow == PL%ncols - 1) field(j, 0)%F = mrb(2 * PL%nx - j + 1)
        field(j,0)%sF = field(j,0)%F

        field(j,0)%Z  = mrb(j + 2 * PL%nx)
        if(PL%myrow == PL%ncols - 1) field(j, 0)%Z = mrb(3 * PL%nx - j + 1)
        field(j,0)%sZ = field(j,0)%Z

      end do
    end if
      
  end do
    
    !print *, "Info sent and recieved in the ", PL%myproc, " processor!"
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    if (PL%myproc /= 0) then
      !print *, "Processor ", PL%myproc, " is up to duty"
      do j = 1, PL%ny
        do i = 1, PL%nx
          index = (j - 1) * PL%nx + i
          soltransf(                    index) = field(i,j)%T
          soltransf(    PL%nx * PL%ny + index) = field(i,j)%F
          soltransf(2 * PL%nx * PL%ny + index) = field(i,j)%Z
        end do
      end do
      ! print *, "Actual size of ", PL%myproc, ":", size(soltransf)
      !print *, "Sending info to main proc"
      !print *, soltransf
      call MPI_SEND(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
      0, 1000 + PL%myproc, MPI_COMM_WORLD, ierr)
      !print *, "Proc" , PL%myproc, "done sending info to main prco"
    else 
      !print *, "Processor 0 is up to duty"
      
      do i = 0, PL%nrows-1
        do j = 0, PL%ncols-1
          
          proc = i * PL%nrows + j
          !print *, "We are in proc", proc, "if u ask proc 0", PL%nrows, PL%ncols
          
          rowsize = i_end_p(proc) - i_str_p(proc) + 1
          colsize = j_end_p(proc) - j_str_p(proc) + 1
          meshsize = rowsize * colsize
          
          if (i == 0 .and. j == 0) then 
            do k = 1,rowsize
              do l = 1,colsize
                
                solution(k ,l, t)%T = field(k,l)%T
                solution(k ,l, t)%F = field(k,l)%F
                solution(k ,l, t)%Z = field(k,l)%Z
              end do
            end do
            
          else
            
            deallocate(soltransf)
            allocate(soltransf( 3 * meshsize))
            ! print *, "Main proc idea of size of proc ", proc, ":", size(soltransf) 
            call MPI_RECV(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
            proc, 1000 + proc, MPI_COMM_WORLD, status, ierr)            
            !print *, "We are chilling in proc 0"
            do l = 1,colsize
              do k = 1,rowsize
                index = (l - 1) * rowsize + k
                
                solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%T = &
                soltransf(               index)
                solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%F = &
                soltransf(    meshsize + index)
                solution(k - 1 + i_str_p(proc),l - 1 + j_str_p(proc), t)%Z = &
                soltransf(2 * meshsize + index)
              end do
            end do
          end if
          
          
        end do
      end do
    end if

    if (PL%myproc == 0 .and. mod(t,param%show) == 0) &
      call write_tstep(solution(:,:,t), param, t)
    
    field0 = field
    if(mod(t,param%show) == 0)  print *, "---------------------------------------------------------------"
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
  end do
  call MPI_FINALIZE(ierr)
end program main