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
  ! use point_m, only: point_t, coefx_1, coefx_2, coefx_3, coefx_4, coefx_5, &
  !                    coefx_6, coefx_7, coefx_8, coefx_9, coefy_1, coefy_2, &
  !                    coefy_3, coefy_4, coefy_5, coefy_6, coefy_7, coefy_8, &
  !                    coefy_9, coefy_10, coefy_11, coefy_12, init_gamma
  use point2_m, only: point_t, init_gamma, wall_exchange_bot_left_boundary, &
                      wall_exchange_bot_boundary, wall_exchange_bot_right_boundary, &
                      wall_exchange_top_left_boundary, wall_exchange_top_boundary, &
                      wall_exchange_top_right_boundary, wall_bot_bd, wall_top_bd
  implicit none

  include 'mpif.h'

  type(param_t) :: param
  type(parll_t) :: PL

  real(DP) :: t_cpu_start, t_cpu_end

  type(point_t), allocatable :: field(:, :), solution(:, :)

  real(DP), allocatable :: d_north_inside_wall(:), d_south_inside_wall(:), d_north_outside_wall(:), d_south_outside_wall(:)
  real(DP) :: error, tmperror
  ! Error for MPI
  integer :: ierr

  ! Integer needed for indeces
  integer :: i, j, m, t, tshow, adi_swp = 1, gauss_seidel_sweps = 1
  integer :: ex_size = 0, iExStart = 0, iExEnd = 0, jSouthOfWall = 0, jNorthOfWall = 0
  !jSouthOfWall, jNorthOfWall = my1, my2 but relative to the processor
  !iExStart, iExEnd, = mx1, mx2 but relative to the processor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialize MPI
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call PL%start_mpi(param, t_cpu_start)

  if (PL%myproc == 0) allocate (solution(param%nx, param%ny))

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
    if (PL%mybot < 0) PL%mybot = MPI_PROC_NULL

    PL%mytop = PL%myproc + PL%ncols
    if (PL%mytop >= PL%nprocs) PL%mytop = MPI_PROC_NULL
  else
    PL%mybot = MPI_PROC_NULL

    PL%mytop = MPI_PROC_NULL
  end if

  if (PL%myproc == 0) call cpu_time(t_cpu_end)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialization of the field
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The 0 index and the n+1 are from processors that are neighbours
  allocate (field(0:PL%nx + 1, 0:PL%ny + 1))

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

      if (wall_bot_bd(field(i, j)%setxc)) then
        ex_size = ex_size + 1

        iExStart = min(iExStart, i)
        if (iExStart == 0) iExStart = i
        ! if (pl%myproc == 1) print *, i, iExStart
        iExEnd = max(iExEnd, i)

        jSouthOfWall = min(j, param%my1)
      else if (wall_top_bd(field(i, j)%setxc)) then
        jNorthOfWall = max(j, param%my2)
      end if
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

  if (iExStart /= 0) then
    do i = iExStart, iExEnd
      do j = jSouthOfWall, jNorthOfWall
        if (wall_bot_bd(field(i, j + 1)%setxc)) then ! my1
          d_north_outside_wall(i) = -(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
                                    /(3.0_DP*param%hy)/(param%bK*param%a**2)!(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
          ! + field(i, j - 2)%T)/(param%hy)/(param%bK*param%a**2)!3.0_DP*((field(i, j)%T + field(i, j + 1)%T)/2.0_DP) - &
          !4.0_DP*((field(i, j - 1)%T + field(i, j)%T)/2.0_DP) + &
          !1.0_DP*((field(i, j - 2)%T + field(i, j - 1)%T)/2.0_DP)/(param%hy*2.0_DP)/(param%bK*param%a**2)

        else if (wall_bot_bd(field(i, j)%setxc)) then ! my1 + 1
          d_south_inside_wall(i) = (13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
                                   /(3.0_DP*param%hy)*(param%bK*param%a**2)!-(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
          !  + field(i, j + 2)%T)*(param%hy)/(param%bK*param%a**2)!-3.0_DP*((field(i, j)%T + field(i, j - 1)%T)/2.0_DP) + &
          !4.0_DP*((field(i, j + 1)%T + field(i, j)%T)/2.0_DP) - &
          !1.0_DP*((field(i, j + 2)%T + field(i, j + 1)%T)/2.0_DP)/(param%hy*2.0_DP)*param%bK*param%a**2

        else if (wall_top_bd(field(i, j)%setxc)) then ! my2-1
          d_north_inside_wall(i) = -(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
                                   /(3.0_DP*param%hy)*(param%bK*param%a**2)!(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
          ! + field(i, j - 2)%T)/(param%hy)*(param%bK*param%a**2)!3.0_DP*((field(i, j)%T + field(i, j + 1)%T)/2.0_DP) - &
          !4.0_DP*((field(i, j - 1)%T + field(i, j)%T)/2.0_DP) + &
          !1.0_DP*((field(i, j - 2)%T + field(i, j - 1)%T)/2.0_DP)/(param%hy*2.0_DP)*param%bK*param%a**2

        else if (wall_top_bd(field(i, j - 1)%setxc)) then ! my2
          d_south_outside_wall(i) = (13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
                                    /(3.0_DP*param%hy)*(param%bK*param%a**2)!-(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
          !  + field(i, j + 2)%T)/(param%hy)/(param%bK*param%a**2)!-3.0_DP*((field(i, j)%T + field(i, j - 1)%T)/2.0_DP) + &
          !4.0_DP*((field(i, j + 1)%T + field(i, j)%T)/2.0_DP) - &
          !1.0_DP*((field(i, j + 2)%T + field(i, j + 1)%T)/2.0_DP)/(param%hy*2.0_DP)/(param%bK*param%a**2)
        end if
        ! if (associated(field(i, j)%setYc, coefy_6)) then ! my1
        !   d_north_outside_wall(i) = (2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
        !                              + field(i, j - 2)%T)/(param%hy)*param%bK/param%a2
        ! else if (associated(field(i, j)%setYc, coefy_12)) then ! my1 + 1
        !   d_south_inside_wall(i) = -(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
        !                              + field(i, j + 2)%T)/(param%hy)/param%bK*param%a2
        ! else if (associated(field(i, j)%setYc, coefy_11)) then ! my2-1
        !   d_north_inside_wall(i) = (2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
        !                             + field(i, j - 2)%T)/(param%hy)/param%bK*param%a2
        ! else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
        !   d_south_outside_wall(i) = -(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
        !                               + field(i, j + 2)%T)/(param%hy)*param%bK/param%a2
        ! end if
      end do
    end do
  end if

  ! CARE: IN CASE THE PROCESSORS END IN MY1 OR MY2 THIS IS A PROBLEM
  ! do i = iExStart, iExEnd
  !   do j = jSouthOfWall, jNorthOfWall
  !     if (associated(field(i, j)%setYc, coefy_6)) then ! my1
  !       field(i, j)%T = 0.5_DP*(field(i, j)%T + field(i, j + 1)%T)
  !       field(i, j + 1)%T = field(i, j)%T
  !     else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
  !       field(i, j)%T = 0.5_DP*(field(i, j)%T + field(i, j - 1)%T)
  !       field(i, j - 1)%T = field(i, j)%T
  !     end if
  !   end do
  ! end do

  call PL%send_info_to_neighbours(field(:, :))
  call PL%send_info_to_proc0(field, solution(:, :))

  if (PL%myproc == 0) call write_tstep(solution(:, :), param, 0)
  tshow = 2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Main Loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do t = 1, param%nt
    m = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Seting coef and system solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, PL%nx
      do j = 1, PL%ny

        call field(i, j)%setXc(param)
        call field(i, j)%setYc(param)

        if (iExStart /= 0) then
        if (wall_bot_bd(field(i, j + 1)%setxc)) then !my1 GAS
          field(i, j)%sT = field(i, j)%sT + d_south_inside_wall(i)*param%hx*field(i, j)%gammaTy
        else if (wall_bot_bd(field(i, j)%setxc)) then ! my1 + 1 SOLID
          field(i, j)%sT = field(i, j)%sT - param%hx*param%a2*d_north_outside_wall(i)
        else if (wall_top_bd(field(i, j - 1)%setxc)) then !my2 GAS
          field(i, j)%sT = field(i, j)%sT - d_north_inside_wall(i)*param%hx*field(i, j)%gammaTy
        else if (wall_top_bd(field(i, j)%setxc)) then ! my2-1 SOLID
          field(i, j)%sT = field(i, j)%sT + param%hx*param%a2*d_south_outside_wall(i)

        end if
        end if
      end do
    end do

    do i = 1, gauss_seidel_sweps
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
    if (iExStart /= 0) then
    do i = iExStart, iExEnd
      do j = jSouthOfWall, jNorthOfWall

        if (wall_bot_bd(field(i, j + 1)%setxc)) then ! my1
          d_north_outside_wall(i) = -(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
                                    /(3.0_DP*param%hy)/(param%bK*param%a**2)!(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
          ! + field(i, j - 2)%T)/(param%hy)/(param%bK*param%a**2)!3.0_DP*((field(i, j)%T + field(i, j + 1)%T)/2.0_DP) - &
          !4.0_DP*((field(i, j - 1)%T + field(i, j)%T)/2.0_DP) + &
          !1.0_DP*((field(i, j - 2)%T + field(i, j - 1)%T)/2.0_DP)/(param%hy*2.0_DP)/(param%bK*param%a**2)

        else if (wall_bot_bd(field(i, j)%setxc)) then ! my1 + 1
          d_south_inside_wall(i) = (13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
                                   /(3.0_DP*param%hy)*(param%bK*param%a**2)!-(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
          !  + field(i, j + 2)%T)*(param%hy)/(param%bK*param%a**2)!-3.0_DP*((field(i, j)%T + field(i, j - 1)%T)/2.0_DP) + &
          !4.0_DP*((field(i, j + 1)%T + field(i, j)%T)/2.0_DP) - &
          !1.0_DP*((field(i, j + 2)%T + field(i, j + 1)%T)/2.0_DP)/(param%hy*2.0_DP)*param%bK*param%a**2

        else if (wall_top_bd(field(i, j)%setxc)) then ! my2-1
          d_north_inside_wall(i) = -(13.0_DP*field(i, j)%T - field(i, j - 1)%T - 4.0_DP*field(i, j + 1)%T) &
                                   /(3.0_DP*param%hy)*(param%bK*param%a**2)!(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
          ! + field(i, j - 2)%T)/(param%hy)*(param%bK*param%a**2)!3.0_DP*((field(i, j)%T + field(i, j + 1)%T)/2.0_DP) - &
          !4.0_DP*((field(i, j - 1)%T + field(i, j)%T)/2.0_DP) + &
          !1.0_DP*((field(i, j - 2)%T + field(i, j - 1)%T)/2.0_DP)/(param%hy*2.0_DP)*param%bK*param%a**2

        else if (wall_top_bd(field(i, j - 1)%setxc)) then ! my2
          d_south_outside_wall(i) = (13.0_DP*field(i, j)%T - field(i, j + 1)%T - 4.0_DP*field(i, j - 1)%T) &
                                    /(3.0_DP*param%hy)/(param%bK*param%a**2)!-(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
          !  + field(i, j + 2)%T)/(param%hy)/(param%bK*param%a**2)!-3.0_DP*((field(i, j)%T + field(i, j - 1)%T)/2.0_DP) + &
          !4.0_DP*((field(i, j + 1)%T + field(i, j)%T)/2.0_DP) - &
          !1.0_DP*((field(i, j + 2)%T + field(i, j + 1)%T)/2.0_DP)/(param%hy*2.0_DP)/(param%bK*param%a**2)
        end if
        ! if (associated(field(i, j)%setYc, coefy_6)) then ! my1
        !   d_north_outside_wall(i) = -(9.0_DP*field(i, j)%T - 8.0_DP*field(i, j + 1)%T &
        !                               - field(i, j - 1)%T)/(3.0_DP*param%hy)*param%bK/param%a2
        ! else if (associated(field(i, j)%setYc, coefy_12)) then ! my1 + 1
        !   d_south_inside_wall(i) = (9.0_DP*field(i, j)%T - 8.0_DP*field(i, j - 1)%T &
        !                             - field(i, j + 1)%T)/(3.0_DP*param%hy)/param%bK*param%a2
        ! else if (associated(field(i, j)%setYc, coefy_11)) then ! my2-1
        !   d_north_inside_wall(i) = -(9.0_DP*field(i, j)%T - 8.0_DP*field(i, j + 1)%T &
        !                              - field(i, j - 1)%T)/(3.0_DP*param%hy)/param%bK*param%a2
        ! else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
        !   d_south_outside_wall(i) = (9.0_DP*field(i, j)%T - 8.0_DP*field(i, j - 1)%T &
        !                              - field(i, j + 1)%T)/(3.0_DP*param%hy)*param%bK/param%a2
        ! end if
        ! if (associated(field(i, j)%setYc, coefy_6)) then ! my1
        !   d_north_outside_wall(i) = (2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
        !                              + field(i, j - 2)%T)/(param%hy)*param%bK/param%a2
        ! else if (associated(field(i, j)%setYc, coefy_12)) then ! my1 + 1
        !   d_south_inside_wall(i) = -(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
        !                              + field(i, j + 2)%T)/(param%hy)/param%bK*param%a2
        ! else if (associated(field(i, j)%setYc, coefy_11)) then ! my2-1
        !   d_north_inside_wall(i) = (2.0_DP*field(i, j)%T - 3.0_DP*field(i, j - 1)%T &
        !                             + field(i, j - 2)%T)/(param%hy)/param%bK*param%a2
        ! else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
        !   d_south_outside_wall(i) = -(2.0_DP*field(i, j)%T - 3.0_DP*field(i, j + 1)%T &
        !                               + field(i, j + 2)%T)/(param%hy)*param%bK/param%a2
        ! end if
      end do
    end do
    end if
    ! CARE: IN CASE THE PROCESSORS END IN MY1 OR MY2 THIS IS A PROBLEM
    ! do i = iExStart, iExEnd
    !   do j = jSouthOfWall, jNorthOfWall
    !     if (associated(field(i, j)%setYc, coefy_6)) then ! my1
    !       field(i, j)%T = 0.5_DP*(field(i, j)%T + field(i, j + 1)%T)
    !       field(i, j + 1)%T = field(i, j)%T
    !     else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
    !       field(i, j)%T = 0.5_DP*(field(i, j)%T + field(i, j - 1)%T)
    !       field(i, j - 1)%T = field(i, j)%T
    !     end if
    !   end do
    ! end do

    ! do i = iExStart, iExEnd
    !   do j = jSouthOfWall, jNorthOfWall
    !     if (associated(field(i, j)%setYc, coefy_6)) then ! my1
    !       field(i, j)%T = 0.5_DP*(field(i, j)%T + field(i, j + 1)%T)
    !       field(i, j + 1)%T = field(i, j)%T
    !     else if (associated(field(i, j)%setYc, coefy_7)) then ! my2
    !       field(i, j)%T = 0.5_DP*(field(i, j)%T + field(i, j - 1)%T)
    !       field(i, j - 1)%T = field(i, j)%T
    !     end if
    !   end do
    ! end do

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

  end do
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)

end program
