program main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Processors numeration
  !  _____________________________________________________________________
  ! |         6         |              7           |           8          |
  ! |___________________|__________________________|______________________|
  ! |         3         |              4           |           5          |
  ! |___________________|__________________________|______________________|
  ! |         0         |              1           |           2          |
  ! |___________________|__________________________|______________________|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use type_m, only: DP
  use param_m, only: param_t
  use paralel_m, only: parll_t
  use coefun_m, only: ADI, write_tstep, gauss_seidel, explicit
  use point_m, only: point_t, coefx_1, coefx_2, coefx_3, coefx_4, coefx_5, &
                     coefx_6, coefx_7, coefx_8, coefx_9, coefy_1, coefy_2, &
                     coefy_3, coefy_4, coefy_5, coefy_6, coefy_7, coefy_8, &
                     coefy_9, coefy_10, coefy_11, coefy_12, init_gamma

  implicit none

  include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'

  type(param_t) :: param
  type(parll_t) :: PL

  real(DP) :: t_cpu_start, t_cpu_end

  type(point_t), allocatable :: field(:, :), field0(:, :), solution(:, :)
  ! Message Send/Recive Left/Riht/Bot/Top
  real(DP), allocatable :: everyerror(:)

  real(DP), allocatable :: d_north_inside_wall(:), d_south_inside_wall(:), d_north_outside_wall(:), d_south_outside_wall(:)
  real(DP) :: error, tmperror
  ! Error for MPI
  integer :: ierr!, ROW_COMM, COL_COMM

  ! Integer needed for indeces
  integer :: i, j, m, t, tshow, ex_size, adi_swp

  adi_swp = 4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MPI SET UP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialize MPI
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call PL%start_mpi(param, t_cpu_start)

  if (PL%myproc == 0) allocate (solution(param%nx, param%ny), &
                                everyerror(PL%nprocs))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Seting who are the neighbours
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Set IDS of neighbours
  PL%myleft = PL%myproc - 1
  if (mod(PL%myproc, PL%ncols) == 0) PL%myleft = MPI_PROC_NULL

  PL%myright = PL%myproc + 1
  if (mod(PL%myright, PL%ncols) == 0) PL%myright = MPI_PROC_NULL

  ! We impose simmetry
  if (PL%nrows /= 1) then
    PL%mybot = PL%myproc - PL%ncols
    if (PL%mybot < 0) PL%mybot = MPI_PROC_NULL !PL%nprocs - PL%myproc - 1

    PL%mytop = PL%myproc + PL%ncols
    if (PL%mytop >= PL%nprocs) PL%mytop = MPI_PROC_NULL !PL%ncols - mod(PL%myproc, PL%ncols) - 1
  else
    PL%mybot = MPI_PROC_NULL

    PL%mytop = MPI_PROC_NULL
  end if

  if (PL%myproc == 0) call cpu_time(t_cpu_end)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialization of the field
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate (field(0:PL%nx + 1, 0:PL%ny + 1))
  allocate (field0(0:PL%nx + 1, 0:PL%ny + 1))

  ex_size = 0
  !print *, "Memory allocated in the ", PL%myproc, " processor!"

  field(:, :)%sT = 0.0_DP
  field(:, :)%sF = 0.0_DP
  field(:, :)%sZ = 0.0_DP

  field(0, 0)%aT(:) = 0.0_DP
  field(0, 0)%aF(:) = 0.0_DP
  field(0, 0)%aZ(:) = 0.0_DP

  field(0, 0)%aT(5) = 1.0_DP
  field(0, 0)%aF(5) = 1.0_DP
  field(0, 0)%aZ(5) = 1.0_DP

  field(PL%nx + 1, 0)%aT(:) = field(0, 0)%aT(:)
  field(PL%nx + 1, 0)%aF(:) = field(0, 0)%aF(:)
  field(PL%nx + 1, 0)%aZ(:) = field(0, 0)%aZ(:)

  do i = 1, PL%nx

    field(i, 0)%aT(:) = field(0, 0)%aT(:)
    field(i, 0)%aF(:) = field(0, 0)%aF(:)
    field(i, 0)%aZ(:) = field(0, 0)%aZ(:)

    do j = 1, PL%ny
      field(0, j)%aT(:) = field(0, 0)%aT(:)
      field(0, j)%aF(:) = field(0, 0)%aF(:)
      field(0, j)%aZ(:) = field(0, 0)%aZ(:)

      field(0, j)%sT = 0.0_DP
      field(0, j)%sF = 1.0_DP
      field(0, j)%sZ = 0.0_DP

      call field(i, j)%set_coef(param, i + PL%istr - 1, &
                                j + PL%jstr - 1)

      call init_gamma(field(i, j), param)
      !if (PL%myproc == 0) print *, i,j, field(i,j)%T, field(i,j)%F, field(i,j)%Z

      if (associated(field(i, j)%setXc, coefy_6)) ex_size = ex_size + 1
      field(PL%nx + 1, j)%aT(:) = field(0, 0)%aT(:)
      field(PL%nx + 1, j)%aF(:) = field(0, 0)%aF(:)
      field(PL%nx + 1, j)%aZ(:) = field(0, 0)%aZ(:)

      field(PL%nx + 1, j)%sT = 0.0_DP
      field(PL%nx + 1, j)%sF = 1.0_DP
      field(PL%nx + 1, j)%sZ = 0.0_DP
    end do

    field(i, PL%ny + 1)%aT(:) = field(0, 0)%aT(:)
    field(i, PL%ny + 1)%aF(:) = field(0, 0)%aF(:)
    field(i, PL%ny + 1)%aZ(:) = field(0, 0)%aZ(:)

  end do

  field(0, PL%ny + 1)%aT(:) = field(0, 0)%aT(:)
  field(PL%nx + 1, PL%ny + 1)%aT(:) = field(0, 0)%aT(:)
  field(0, PL%ny + 1)%aF(:) = field(0, 0)%aF(:)
  field(PL%nx + 1, PL%ny + 1)%aF(:) = field(0, 0)%aF(:)
  field(0, PL%ny + 1)%aZ(:) = field(0, 0)%aZ(:)
  field(PL%nx + 1, PL%ny + 1)%aZ(:) = field(0, 0)%aZ(:)

  allocate (d_north_inside_wall(param%nx), d_south_inside_wall(param%nx), &
            d_north_outside_wall(param%nx), d_south_outside_wall(param%nx))
  do i = 1, PL%nx
    do j = 1, PL%ny
      if (associated(field(i, j)%setYc, coefy_6)) then ! my1
        d_north_outside_wall(i) = -(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
                                  /(3.0_DP*param%hy)*param%bK
      else if (associated(field(i, j)%setYc, coefy_12)) then ! my1 + 1
        d_south_inside_wall(i) = (13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
                                 /(3.0_DP*param%hy)/param%bK
      else if (associated(field(i, j)%setYc, coefy_11)) then ! my2-1
        d_north_inside_wall(i) = -(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
                                 /(3.0_DP*param%hy)/param%bK
      else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
        d_south_outside_wall(i) = (13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
                                  /(3.0_DP*param%hy)*param%bK
      end if
    end do
  end do
  field0 = field

  call PL%send_info_to_neighbours(field(:, :))
  call PL%send_info_to_proc0(field, solution(:, :))

  if (PL%myproc == 0) call write_tstep(solution(:, :), param, 0)

  tshow = 2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Main Loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do t = 1, param%nt
    m = 0
    !print *, t
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Seting coef and system solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print *, "Processor:", PL%myproc, param%nt, PL%nx, PL%ny, param%my1, param%my2
    do i = 1, PL%nx
      do j = 1, PL%ny
        if (.not. associated(field(i, j)%setXc)) print *, i, j
        call field(i, j)%setXc(param)
        call field(i, j)%setYc(param)
        ! if (i == 69 .and. j == 45) print *, "-", field(69,44)%T

        if (associated(field(i, j)%setYc, coefy_6)) then !my1
          field(i, j)%sT = field(i, j)%sT + d_south_inside_wall(i)*param%hx*field(i, j)%gammaTy
        else if (associated(field(i, j)%setYc, coefy_12)) then ! my1 + 1
          field(i, j)%sT = field(i, j)%sT - d_north_outside_wall(i)*param%hx*param%alpha*param%a2
        else if (associated(field(i, j)%setYc, coefy_7)) then !my2
          field(i, j)%sT = field(i, j)%sT - d_north_inside_wall(i)*param%hx*field(i, j)%gammaTy
        else if (associated(field(i, j)%setYc, coefy_11)) then ! my2-1
          field(i, j)%sT = field(i, j)%sT + d_south_outside_wall(i)*param%hx*param%alpha*param%a2

        end if
      end do
    end do

    !print *, "Coeficien functions done in the ", PL%myproc, " processor!"
    call explicit(field(0:PL%nx + 1, 0:PL%ny + 1)%aT(1), field(0:PL%nx + 1, 0:PL%ny + 1)%aT(2), &
                  field(0:PL%nx + 1, 0:PL%ny + 1)%aT(3), field(0:PL%nx + 1, 0:PL%ny + 1)%aT(4), &
                  field(0:PL%nx + 1, 0:PL%ny + 1)%aT(5), field(:, :)%sT, field0(:, :)%T, &
                  field(:, :)%T)
    field%T = field0%T + field%T*param%ht/(param%hx*param%hy)
    error = norm2(field(:, :)%T - field0(:, :)%T)

    call explicit(field(0:PL%nx + 1, 0:PL%ny + 1)%aF(1), field(0:PL%nx + 1, 0:PL%ny + 1)%aF(2), &
                  field(0:PL%nx + 1, 0:PL%ny + 1)%aF(3), field(0:PL%nx + 1, 0:PL%ny + 1)%aF(4), &
                  field(0:PL%nx + 1, 0:PL%ny + 1)%aF(5), field(:, :)%sF, field0(:, :)%F, &
                  field(:, :)%F)
    field%F = field0%F + field%F*param%ht/(param%hx*param%hy)
    error = max(norm2(field(:, :)%F - field0(:, :)%F), error)

    call explicit(field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(1), field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(2), &
                  field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(3), field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(4), &
                  field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(5), field(:, :)%sZ, field0(:, :)%Z, &
                  field(:, :)%Z)
    field%Z = field0%Z + field%Z*param%ht/(param%hx*param%hy)
    error = max(norm2(field(:, :)%Z - field0(:, :)%Z), error)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Setting exchange boundaries
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, PL%nx
      do j = 1, PL%ny
        if (associated(field(i, j)%setYc, coefy_6)) then ! my1
          d_north_outside_wall(i) = (4.0_DP*field(i, j)%T - 6.0_DP*field(i, j - 1)%T &
                                     + 2.0_DP*field(i, j - 2)%T)/(3.0_DP*param%hy)*param%bK!-(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
          !/(3.0_DP*param%hy)*param%bK
        else if (associated(field(i, j)%setYc, coefy_12)) then ! my1 + 1
          d_south_inside_wall(i) = -(4.0_DP*field(i, j)%T - 6.0_DP*field(i, j + 1)%T &
                                     + 2.0_DP*field(i, j + 2)%T)/(3.0_DP*param%hy)/param%bK!(13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
          !/(3.0_DP*param%hy)/param%bK
        else if (associated(field(i, j)%setYc, coefy_11)) then ! my2-1
          d_north_inside_wall(i) = (4.0_DP*field(i, j)%T - 6.0_DP*field(i, j - 1)%T &
                                    + 2.0_DP*field(i, j - 2)%T)/(3.0_DP*param%hy)/param%bK!-(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
          !/(3.0_DP*param%hy)/param%bK
        else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
          d_south_outside_wall(i) = -(4.0_DP*field(i, j)%T - 6.0_DP*field(i, j + 1)%T &
                                      + 2.0_DP*field(i, j + 2)%T)/(3.0_DP*param%hy)*param%bK!(13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
          !/(3.0_DP*param%hy)*param%bK
        end if
      end do
    end do
    ! open(555, file = "help.txt")
    ! do i = 1, param%nx
    !   write(555, *) ex_wbot_bd(i), ex_wtop_bd(param%nx-i+1)

    ! end do
    ! close(555)
    ! stop
    call PL%send_info_to_neighbours(field(:, :))

    call PL%send_errors(error)
    if (mod(t, param%show) == 0) then

      if (PL%myproc == 0) then
        call PL%send_info_to_proc0(field(:, :), solution(:, :))
        call write_tstep(solution(:, :), param, tshow)
        print *, "---------------------------------------------------------------"
        print *, "M = ", t, error
        print *, "---------------------------------------------------------------"
      else
        call PL%send_info_to_proc0(field(:, :))
      end if
      tshow = int(t/param%show) + 1
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    field0 = field
  end do
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)

contains

  subroutine stopit(a)
    real(DP) :: a
    if (a /= 0) then
      print *, "here"
      stop
    else
      print *, "not here"
    end if
  end subroutine stopit
end program
