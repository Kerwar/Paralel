module paralel_m

  use type_m, only: DP
  use param_m, only: param_t
  use point_m, only: point_t

  implicit none

  include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'

  type, public :: parll_t
    ! Number of processor and nuber of this processor
    integer :: nprocs, myproc
    ! Number of columns and rows of the decomposition
    integer :: nrows, ncols
    ! Number of the column and row of this processor
    integer :: myrow, mycol
    ! Indeces where this processor has its calculations
    integer :: istr, iend, jstr, jend
    ! Number of processors that have an extra point
    integer :: nxadd, nyadd
    ! IDS for row and col COMMS and size of them
    integer :: myid_row, myid_col, nodes_row, nodes_col
    ! IDS of left and right processors, top and bot
    integer :: myleft, myright, mytop, mybot
    ! Size of the mesh in this processor
    integer :: nx, ny

    ! Indeces for data treatment in proc0
    integer, allocatable :: i_str_p(:), i_end_p(:), &
                            j_str_p(:), j_end_p(:)
  contains

    procedure :: start_mpi
    procedure :: send_info_to_neighbours
    procedure :: send_info_to_proc0
    procedure :: send_errors

    procedure, private :: test
  end type parll_t

contains

  subroutine start_mpi(this, param, t_cpu_start)
    class(parll_t) :: this
    type(param_t)  :: param
    real(DP)       :: t_cpu_start
    integer        :: ierr
    integer :: i, j

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, this%nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, this%myproc, ierr)

    if (this%myproc == 0) call cpu_time(t_cpu_start)

    call param%read_param()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Domain decomposition in processors
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    this%nrows = int(sqrt(float(param%ny)/float(param%nx)* &
                          float(this%nprocs)))
    if (this%nrows == 0) this%nrows = 1

    this%ncols = int(float(this%nprocs)/float(this%nrows))

    do while (this%nrows*this%ncols /= this%nprocs)
      this%nrows = this%nrows - 1
      this%ncols = int(float(this%nprocs)/float(this%nrows))
    end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Domain decomposition and setting meshing sizes
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !NEED TO CHECK IF THIS DECOMPOSITION WORKS FOR MORE THAN 2 PROC
    this%nxadd = mod(param%nx, this%ncols)
    this%nyadd = mod(param%ny, this%nrows)

    this%myrow = min(int(float(this%myproc)/float(this%ncols)), this%ncols - 1)
    this%mycol = this%myproc - this%myrow*this%ncols

    if (this%mycol < this%nxadd) then
      this%istr = this%mycol* &
                  (int(float(param%nx)/float(this%ncols)) + 1) + 1
      this%iend = (this%mycol + 1)* &
                  (int(float(param%nx)/float(this%ncols)) + 1)
    else
      this%istr = this%nxadd + this%mycol* &
                  int(float(param%nx)/float(this%ncols)) + 1
      this%iend = this%nxadd + (this%mycol + 1)* &
                  int(float(param%nx)/float(this%ncols))
    end if

    if (this%myrow < this%nyadd) then
      this%jstr = this%myrow* &
                  (int(float(param%ny)/float(this%nrows)) + 1) + 1
      this%jend = (this%myrow + 1)* &
                  (int(float(param%ny)/float(this%nrows)) + 1)
    else
      this%jstr = this%nyadd + this%myrow* &
                  int(float(param%ny)/float(this%nrows)) + 1
      this%jend = this%nyadd + (this%myrow + 1)* &
                  int(float(param%ny)/float(this%nrows))
    end if
    call this%test()
    this%nx = this%iend - this%istr + 1
    this%ny = this%jend - this%jstr + 1

    if (this%myproc == 0) then
      allocate (this%i_str_p(0:this%nprocs - 1), this%i_end_p(0:this%nprocs - 1), &
                this%j_str_p(0:this%nprocs - 1), this%j_end_p(0:this%nprocs - 1))

      do j = 0, this%nrows - 1
        do i = 0, this%nxadd - 1
          this%i_str_p(i + j*this%ncols) = i*(param%nx/this%ncols + 1) + 1
          this%i_end_p(i + j*this%ncols) = (i + 1)*(param%nx/this%ncols + 1)
        enddo
        do i = this%nxadd, this%ncols - 1
          this%i_str_p(i + j*this%ncols) = this%nxadd + &
                                           i*(param%nx/this%ncols) + 1
          this%i_end_p(i + j*this%ncols) = this%nxadd + &
                                           (i + 1)*(param%nx/this%ncols)
        enddo
      enddo

      do i = 0, this%ncols - 1
        do j = 0, this%nyadd - 1
          this%j_str_p(j*this%ncols + i) = j*(param%ny/this%nrows + 1) + 1
          this%j_end_p(j*this%ncols + i) = (j + 1)*(param%ny/this%nrows + 1)
        enddo
        do j = this%nyadd, this%nrows - 1
          this%j_str_p(j*this%ncols + i) = this%nyadd + &
                                           j*(param%ny/this%nrows) + 1
          this%j_end_p(j*this%ncols + i) = this%nyadd + &
                                           (j + 1)*(param%ny/this%nrows)
        enddo
      enddo
    end if
  end subroutine start_mpi

  subroutine send_info_to_neighbours(this, field)
    class(parll_t), intent(in) :: this
    type(point_t), intent(inout) :: field(0:, 0:)

    real(DP), allocatable :: msl(:), mrl(:), msr(:), mrr(:), &
                             msb(:), mrb(:), mst(:), mrt(:)
    integer :: i, j

    integer, parameter :: LEFT_TAG = 10, RIGHT_TAG = 12, TOP_TAG = 13, BOT_TAG = 14
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)

    allocate (msl(3*this%ny), mrl(3*this%ny))
    allocate (msr(3*this%ny), mrr(3*this%ny))
    allocate (msb(3*this%nx), mrb(3*this%nx))
    allocate (mst(3*this%nx), mrt(3*this%nx))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Sending info left
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j = 1, this%ny
      msl(j) = field(1, j)%T
      msl(j + this%ny) = field(1, j)%F
      msl(j + 2*this%ny) = field(1, j)%Z
    end do

    call MPI_SENDRECV(msl, size(msl), MPI_DOUBLE_PRECISION, &
                      this%myleft, LEFT_TAG, mrr, size(mrr), MPI_DOUBLE_PRECISION, &
                      this%myright, LEFT_TAG, MPI_COMM_WORLD, status, ierr)
    do j = 1, this%ny
      field(this%nx + 1, j)%T = mrr(j)
      field(this%nx + 1, j)%sT = field(this%nx + 1, j)%T

      field(this%nx + 1, j)%F = mrr(j + this%ny)
      field(this%nx + 1, j)%sF = field(this%nx + 1, j)%F

      field(this%nx + 1, j)%Z = mrr(j + 2*this%ny)
      field(this%nx + 1, j)%sZ = field(this%nx + 1, j)%Z

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Sending info right
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j = 1, this%ny
      msr(j) = field(this%nx, j)%T
      msr(j + this%ny) = field(this%nx, j)%F
      msr(j + 2*this%ny) = field(this%nx, j)%Z
    end do

    call MPI_SENDRECV(msr, size(msr), MPI_DOUBLE_PRECISION, &
                      this%myright, RIGHT_TAG, mrl, size(mrl), MPI_DOUBLE_PRECISION, &
                      this%myleft, RIGHT_TAG, MPI_COMM_WORLD, status, ierr)

    do j = 1, this%ny
      field(0, j)%T = mrl(j)
      field(0, j)%sT = field(0, j)%T

      field(0, j)%F = mrl(j + this%ny)
      field(0, j)%sF = field(0, j)%F

      field(0, j)%Z = mrl(j + 2*this%ny)
      field(0, j)%sZ = field(0, j)%Z
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Sending info bot
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j = 1, this%nx
      msb(j) = field(j, 1)%T
      msb(j + this%nx) = field(j, 1)%F
      msb(j + 2*this%nx) = field(j, 1)%Z
    end do

    call MPI_SENDRECV(msb, size(msb), MPI_DOUBLE_PRECISION, &
                      this%mybot, BOT_TAG, mrt, size(mrt), MPI_DOUBLE_PRECISION, &
                      this%mytop, BOT_TAG, MPI_COMM_WORLD, status, ierr)

    do j = 1, this%nx
      field(j, this%ny + 1)%T = mrt(j)
      if (this%mytop == MPI_PROC_NULL) field(j, this%ny + 1)%T = field(j, this%ny)%T !mrt(this%nx - j + 1)
      field(j, this%ny + 1)%sT = field(j, this%ny + 1)%T

      field(j, this%ny + 1)%F = mrt(j + this%nx)
      if (this%mytop == MPI_PROC_NULL) field(j, this%ny + 1)%F = field(j, this%ny)%F !mrt(2 * this%nx - j + 1)
      field(j, this%ny + 1)%sF = field(j, this%ny + 1)%F

      field(j, this%ny + 1)%Z = mrt(j + 2*this%nx)
      if (this%mytop == MPI_PROC_NULL) field(j, this%ny + 1)%Z = field(j, this%ny)%Z !mrt(3 * this%nx - j + 1)
      field(j, this%ny + 1)%sZ = field(j, this%ny + 1)%Z

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Sending info top
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j = 1, this%nx
      mst(j) = field(j, this%ny)%T
      mst(j + this%nx) = field(j, this%ny)%F
      mst(j + 2*this%nx) = field(j, this%ny)%Z
    end do

    call MPI_SENDRECV(mst, size(mst), MPI_DOUBLE_PRECISION, &
                      this%mytop, TOP_TAG, mrb, size(mrb), MPI_DOUBLE_PRECISION, &
                      this%mybot, TOP_TAG, MPI_COMM_WORLD, status, ierr)
    do j = 1, this%nx
      field(j, 0)%T = mrb(j)
      if (this%mybot == MPI_PROC_NULL) field(j, 0)%T = field(j, 1)%T !mrb(this%nx - j + 1)
      field(j, 0)%sT = field(j, 0)%T

      field(j, 0)%F = mrb(j + this%nx)
      if (this%mybot == MPI_PROC_NULL) field(j, 0)%F = field(j, 1)%F !mrb(2 * this%nx - j + 1)
      field(j, 0)%sF = field(j, 0)%F

      field(j, 0)%Z = mrb(j + 2*this%nx)
      if (this%mybot == MPI_PROC_NULL) field(j, 0)%Z = field(j, 1)%Z !mrb(3 * this%nx - j + 1)
      field(j, 0)%sZ = field(j, 0)%Z

    end do

  end subroutine send_info_to_neighbours

  subroutine send_info_to_proc0(this, field, solution)
    class(parll_t), intent(in) :: this
    type(point_t), intent(in) :: field(0:, 0:)
    type(point_t), intent(inout), optional :: solution(:, :)

    real(DP), allocatable :: soltransf(:)
    integer :: i, j, k, l
    integer :: meshsize, rowsize, colsize, index, proc
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Sending info to proc 0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (this%myproc /= 0) then

      allocate (soltransf(3*this%nx*this%ny))
      do j = 1, this%ny
        do i = 1, this%nx
          index = (j - 1)*this%nx + i
          soltransf(index) = field(i, j)%T
          soltransf(this%nx*this%ny + index) = field(i, j)%F
          soltransf(2*this%nx*this%ny + index) = field(i, j)%Z
        end do
      end do

      call MPI_SEND(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
                    0, 1000 + this%myproc, MPI_COMM_WORLD, ierr)

      deallocate (soltransf)
    else

      do i = 0, this%nrows - 1
        do j = 0, this%ncols - 1

          proc = i*this%nrows + j

          rowsize = this%i_end_p(proc) - this%i_str_p(proc) + 1
          colsize = this%j_end_p(proc) - this%j_str_p(proc) + 1
          meshsize = rowsize*colsize

          if (i == 0 .and. j == 0) then

            do k = 1, rowsize
              do l = 1, colsize
                solution(k, l)%T = field(k, l)%T
                solution(k, l)%F = field(k, l)%F
                solution(k, l)%Z = field(k, l)%Z
              end do
            end do
          else

            allocate (soltransf(3*meshsize))

            call MPI_RECV(soltransf, size(soltransf), MPI_DOUBLE_PRECISION, &
                          proc, 1000 + proc, MPI_COMM_WORLD, status, ierr)

            do l = 1, colsize
              do k = 1, rowsize
                index = (l - 1)*rowsize + k
                solution(k - 1 + this%i_str_p(proc), l - 1 + this%j_str_p(proc))%T = &
                  soltransf(index)
                solution(k - 1 + this%i_str_p(proc), l - 1 + this%j_str_p(proc))%F = &
                  soltransf(meshsize + index)
                solution(k - 1 + this%i_str_p(proc), l - 1 + this%j_str_p(proc))%Z = &
                  soltransf(2*meshsize + index)
              end do
            end do

            deallocate (soltransf)
          end if

        end do
      end do
    end if
  end subroutine send_info_to_proc0

  subroutine send_errors(this, error)
    class(parll_t), intent(in) :: this
    real(DP), intent(inout) :: error

    integer :: i, j
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    real(DP) :: proc, tmperror

    if (this%myproc /= 0) then
      call MPI_SEND(error, 1, MPI_DOUBLE_PRECISION, &
                    0, 2000 + this%myproc, MPI_COMM_WORLD, ierr)
    else if (this%myproc == 0) then

      do i = 0, this%nrows - 1
        do j = 0, this%ncols - 1

          proc = i*this%nrows + j
          if (i /= 0 .and. j /= 0) then
            call MPI_RECV(tmperror, 1, MPI_DOUBLE_PRECISION, &
                          proc, 2000 + proc, MPI_COMM_WORLD, status, ierr)

            error = max(error, tmperror)
          end if

        end do
      end do
    end if

  end subroutine send_errors

  subroutine test(this)
    class(parll_t) :: this

    print "(A16, I2, A17, I2, A11)", "My processor is:", this%myproc, " from the pool of", &
      this%nprocs, " processors"
  end subroutine test
end module paralel_m
