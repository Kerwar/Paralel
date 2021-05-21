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
  use coefun_m, only: ADI, write_tstep, gauss_seidel
  use point_m, only: point_t, coefx_1, coefx_2, coefx_3, coefx_4, coefx_5, &
                     coefx_6, coefx_7, coefx_8, coefx_9, coefy_1, coefy_2, &
                     coefy_3, coefy_4, coefy_5, coefy_6, coefy_7, coefy_8, &
                     coefy_9, coefy_10, coefy_11, coefy_12, init_gamma

  implicit none

  include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'

  type(param_t) :: param
  type(parll_t) :: PL

  real(DP) :: t_cpu_start, t_cpu_end

  type(point_t), allocatable :: field(:, :), field0(:, :), solution(:, :, :)
  ! Message Send/Recive Left/Riht/Bot/Top
  real(DP), allocatable :: everyerror(:)

  real(DP), allocatable :: ex_wtop_bd(:), ex_wbot_bd(:), ex_top_bd(:), ex_bot_bd(:)
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

  if (PL%myproc == 0) allocate (solution(param%nx, param%ny, int(param%nt/param%show) + 1), &
                                everyerror(PL%nprocs))

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

  allocate (ex_wtop_bd(param%nx), ex_wbot_bd(param%nx), ex_top_bd(param%nx), ex_bot_bd(param%nx))
  do i = 1, PL%nx
    do j = 1, PL%ny
      if (associated(field(i, j)%setYc, coefy_6)) then
        ex_wbot_bd(i) = (9.0_DP*field(i, j)%T - field(i, j - 1)%T - 8.0_DP*field(i, j + 1)%T) &
                        /(3.0_DP*param%hy)/param%bK
        ex_top_bd(i) = -(9.0_DP*field(i, j + 1)%T - field(i, j + 2)%T - 8.0_DP*field(i, j)%T) &
                       /(3.0_DP*param%hy)*param%bK
      else if (associated(field(i, j)%setYc, coefy_7)) then
        ex_wtop_bd(i) = -(9.0_DP*field(i, j)%T - field(i, j + 1)%T - 8.0_DP*field(i, j - 1)%T) &
                        /(3.0_DP*param%hy)/param%bK
        ex_bot_bd(i) = (9.0_DP*field(i, j - 1)%T - field(i, j - 2)%T - 8.0_DP*field(i, j)%T) &
                       /(3.0_DP*param%hy)*param%bK
      end if
    end do
  end do
  field0 = field

  ! if (pl%myproc == 1) then
  !   open (333, file="coef")
  !   do i = 1, pl%nx
  !     if (associated(field(i, param%my2 - 1)%setyc, coefy_11)) then
  !       write (333, *) "Yes"
  !     else if (associated(field(i, param%my2 - 1)%setyc, coefy_5)) then
  !     else
  !       write (333, *) i
  !     end if
  !   end do
  !   close (333)
  ! end if
  call PL%send_info_to_neighbours(field(:, :))
  call PL%send_info_to_proc0(field, solution(:, :, 1))

  if (PL%myproc == 0) call write_tstep(solution(:, :, 1), param, 0)
  call stopit(field(101, 86)%sT)
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
        call stopit(field(101, 86)%sT)
        call field(i, j)%setYc(param)

        call stopit(field(101, 86)%sT)
        field(i, j)%sT = field(i, j)%sT + field0(i, j)%T*param%txy
        field(i, j)%sF = field(i, j)%sF + field0(i, j)%F*param%txy
        field(i, j)%sZ = field(i, j)%sZ + field0(i, j)%Z*param%txy

        if (associated(field(i, j)%setYc, coefy_6)) then
          !field(i, j + 1)%sT = field(i, j)%sT - ex_wbot_bd(i)*param%hx*param%alpha*param%a2
          !field(i, j)%sT = field(i, j)%sT + ex_top_bd(i)*param%hx*field(i, j)%gammaTy
        else if (associated(field(i, j)%setYc, coefy_7)) then
          !field(i, j)%sT = field(i, j)%sT - ex_bot_bd(i)*param%hx*field(i, j)%gammaTy
          !field(i, j - 1)%sT = field(i, j)%sT + ex_wtop_bd(i)*param%hx*param%alpha*param%a2
        else if (associated(field(i, j)%setYc, coefy_8)) then
          field(i, j)%sT = field(i, j)%sT + param%m*param%hx*field(i, j)%v(4)*field(i, j + 1)%T
          field(i, j)%sF = field(i, j)%sF + param%m*param%hx*field(i, j)%v(4)*field(i, j + 1)%F
          field(i, j)%sZ = field(i, j)%sZ + param%m*param%hx*field(i, j)%v(4)*field(i, j + 1)%Z
        else if (associated(field(i, j)%setYc, coefy_9)) then
          field(i, j)%sT = field(i, j)%sT + param%m*param%hx*field(i, j)%v(3)*field(i, j - 1)%T
          field(i, j)%sF = field(i, j)%sF + param%m*param%hx*field(i, j)%v(3)*field(i, j - 1)%F
          field(i, j)%sZ = field(i, j)%sZ + param%m*param%hx*field(i, j)%v(3)*field(i, j - 1)%Z
        end if
      end do
    end do
    open (555, file="help.txt")
    do i = param%mx1, param%mx2
    do j = param%my1 + 1, param%my2 - 1
      if (associated(field(i, j)%setyc, coefy_4)) then
        write (555, "(2I4, 5F14.8, A1, F14.8)") i, j, field(i, j)%aT(:), "|", field(i, j)%sT
      else if (associated(field(i, j)%setyc, coefy_12)) then
        write (555, "(2I4, 5F14.8, A1, F14.8, A1)") i, j, field(i, j)%aT(:), "|", field(i, j)%sT, "-"
      else if (associated(field(i, j)%setyc, coefy_11)) then
        write (555, "(2I4, 5F14.8, A1, F14.8, A1)") i, j, field(i, j)%aT(:), "|", field(i, j)%sT, "+"
      else
        write (555, *) "Nany dafuck"
      end if
    end do
    end do
    close (555)
    if (t == 2) stop
    !if (PL%myproc == 1 ) then

    !end if
    !print *, "Coeficien functions done in the ", PL%myproc, " processor!"
    do i = 1, 5
      call gauss_seidel(field(0:PL%nx + 1, 0:PL%ny + 1)%aT(1), field(0:PL%nx + 1, 0:PL%ny + 1)%aT(2), &
                        field(0:PL%nx + 1, 0:PL%ny + 1)%aT(3), field(0:PL%nx + 1, 0:PL%ny + 1)%aT(4), &
                        field(0:PL%nx + 1, 0:PL%ny + 1)%aT(5), field(:, :)%sT, field(:, :)%T, &
                        tmperror, adi_swp, param%valpha)
      error = max(tmperror, 0.0_DP)

      call gauss_seidel(field(0:PL%nx + 1, 0:PL%ny + 1)%aF(1), field(0:PL%nx + 1, 0:PL%ny + 1)%aF(2), &
                        field(0:PL%nx + 1, 0:PL%ny + 1)%aF(3), field(0:PL%nx + 1, 0:PL%ny + 1)%aF(4), &
                        field(0:PL%nx + 1, 0:PL%ny + 1)%aF(5), field(:, :)%sF, field(:, :)%F, &
                        tmperror, adi_swp, param%valpha)
      error = max(tmperror, error)

      call gauss_seidel(field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(1), field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(2), &
                        field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(3), field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(4), &
                        field(0:PL%nx + 1, 0:PL%ny + 1)%aZ(5), field(:, :)%sZ, field(:, :)%Z, &
                        tmperror, adi_swp, param%valpha)
      error = max(tmperror, error)

      call PL%send_info_to_neighbours(field(:, :))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Setting exchange boundaries
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, PL%nx
      do j = 1, PL%ny
        if (associated(field(i, j)%setYc, coefy_6)) then
          ex_wbot_bd(i) = -(9.0_DP*field(i, j)%T - field(i, j - 1)%T - 8.0_DP*field(i, j + 1)%T) &
                          /(3.0_DP*param%hy)/param%bK
          ex_top_bd(i) = (9.0_DP*field(i, j + 1)%T - field(i, j + 2)%T - 8.0_DP*field(i, j)%T) &
                         /(3.0_DP*param%hy)*param%bK
        else if (associated(field(i, j)%setYc, coefy_7)) then
          ex_wtop_bd(i) = (9.0_DP*field(i, j)%T - field(i, j + 1)%T - 8.0_DP*field(i, j - 1)%T) &
                          /(3.0_DP*param%hy)/param%bK
          ex_bot_bd(i) = -(9.0_DP*field(i, j - 1)%T - field(i, j - 2)%T - 8.0_DP*field(i, j)%T) &
                         /(3.0_DP*param%hy)*param%bK
        end if
      end do
    end do
    call PL%send_info_to_neighbours(field(:, :))

    call PL%send_errors(error)
    if (mod(t, param%show) == 0) then

      if (PL%myproc == 0) then
        call PL%send_info_to_proc0(field(:, :), solution(:, :, tshow))
        call write_tstep(solution(:, :, tshow), param, tshow)
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
    end if
  end subroutine stopit
end program
