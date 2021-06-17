module point2_m

  use type_m, only: DP
  use param_m, only: param_t
  use interfaces_m, only: coef

  implicit none

  type, public :: point_t
    real(DP) :: T, F, Z
    ! u_w, u_e, u_s, u_n
    real(DP) :: v(4)
    ! W, E, S, N, P
    real(DP) :: aT(5), aF(5), aZ(5), sT, sF, sZ

    real(DP) :: gammaTx, gammaTy
    real(DP) :: gammaFx, gammaFy
    real(DP) :: gammaZx, gammaZy
    procedure(coef), nopass, pointer :: cT, cF, cZ
    procedure(co), pointer :: setXc, setYc

  contains

    procedure :: set_coef
    !procedure :: gamma => init_gamma
  end type point_t

  abstract interface
    subroutine co(this, param)
      import :: DP, param_t, point_t
      class(point_t), intent(inout) :: this
      type(param_t), intent(in) :: param
    end subroutine
  end interface

contains

  subroutine init_gamma(this, param)
    type(point_t), intent(inout) :: this
    type(param_t), intent(in) :: param

    this%gammaTx = 1.0_DP

    this%gammaTy = param%a2

    this%gammaFx = 1.0_DP/param%LeF

    this%gammaFy = param%a2/param%LeF

    this%gammaZx = 1.0_DP/param%LeZ

    this%gammaZy = param%a2/param%LeZ

  end subroutine init_gamma

  ! --------------------------------------------
  !                    ---->
  ! -------------------------------------------- my2 does not belong to the wall
  !                --------------                mx1 and mx2 belong to the wall
  ! -------------------------------------------- my1 does not belongs to the wall
  !                    <----
  ! --------------------------------------------

  subroutine set_coef(this, param, nx, ny)
    class(point_t), intent(inout) :: this
    type(param_t), intent(in) :: param
    integer, intent(in) :: nx, ny
    real(DP)       :: rr, pu, pd

    if (nx == 1) then
      if (ny == 1) then

        this%setXc => right_flow_channel_bot_simmetry_outlet
        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > 1 .and. ny < param%my1) then
        this%setXc => right_flow_channel_interior
        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1) then
        this%setXc => right_flow_channel_top_wall_outlet

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my2) then
        this%setXc => left_flow_channel_bot_wall_inlet

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > param%my2 .and. ny < param%ny) then
        this%setXc => left_flow_channel_inlet

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%ny) then
        this%setXc => left_flow_channel_top_simmetry_inlet

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else
        this%setXc => wall_no_exchange

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      end if

    else if (nx == param%nx) then
      if (ny == 1) then
        this%setXc => right_flow_channel_bot_simmetry_inlet

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > 1 .and. ny < param%my1) then
        this%setXc => right_flow_channel_inlet

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1) then
        this%setXc => right_flow_channel_top_wall_inlet

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my2) then
        this%setXc => left_flow_channel_bot_wall_outlet

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > param%my2 .and. ny < param%ny) then
        this%setXc => left_flow_channel_outlet

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%ny) then
        this%setXc => left_flow_channel_top_simmetry_outlet

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else
        this%setXc => wall_no_exchange

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      end if

    else if (nx == param%mx1) then
      if (ny == 1) then
        this%setXc => right_flow_channel_bot_simmetry

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > 1 .and. ny < param%my1) then
        this%setXc => right_flow_channel_interior

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1) then
        this%setXc => right_flow_channel_top_wall

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my2) then
        this%setXc => left_flow_channel_bot_wall

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > param%my2 .and. ny < param%ny) then
        this%setXc => left_flow_channel_interior

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%ny) then
        this%setXc => left_flow_channel_top_simmetry

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1 + 1) then
        this%setXc => wall_exchange_bot_left_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      else if (ny == param%my2 - 1) then
        this%setXc => wall_exchange_top_left_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      else
        this%setXc => wall_exchange_left_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      end if

    else if (nx == param%mx2) then
      if (ny == 1) then
        this%setXc => right_flow_channel_bot_simmetry

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > 1 .and. ny < param%my1) then
        this%setXc => right_flow_channel_interior

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1) then
        this%setXc => right_flow_channel_top_wall

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my2) then
        this%setXc => left_flow_channel_bot_wall

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > param%my2 .and. ny < param%ny) then
        this%setXc => left_flow_channel_interior

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%ny) then
        this%setXc => left_flow_channel_top_simmetry

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1 + 1) then
        this%setXc => wall_exchange_bot_right_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      else if (ny == param%my2 - 1) then
        this%setXc => wall_exchange_top_right_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      else
        this%setXc => wall_exchange_right_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      end if

    else if (nx > param%mx1 .and. nx < param%mx2) then
      if (ny == 1) then
        this%setXc => right_flow_channel_bot_simmetry

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > 1 .and. ny < param%my1) then
        this%setXc => right_flow_channel_interior

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1) then
        this%setXc => right_flow_channel_top_wall

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my2) then
        this%setXc => left_flow_channel_bot_wall

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > param%my2 .and. ny < param%ny) then
        this%setXc => left_flow_channel_interior

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%ny) then
        this%setXc => left_flow_channel_top_simmetry

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1 + 1) then
        this%setXc => wall_exchange_bot_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      else if (ny == param%my2 - 1) then
        this%setXc => wall_exchange_top_boundary

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      else
        this%setXc => wall_exchange

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      end if

    else
      if (ny == 1) then
        this%setXc => right_flow_channel_bot_simmetry

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > 1 .and. ny < param%my1) then
        this%setXc => right_flow_channel_interior

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my1) then
        this%setXc => right_flow_channel_top_wall

        pu = param%hy*(-ny + 0.5_DP + param%my1)
        pd = param%hy*(2.0_DP*param%my1 + 1.0_DP) - pu
        this%v(1) = -6.0_DP*param%m*pu*pd
        this%v(2) = -6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax + param%xhs)**2 + &
                  ((ny - 1)*param%hy - param%yhs)**2)

        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%my2) then
        this%setXc => left_flow_channel_bot_wall

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny > param%my2 .and. ny < param%ny) then
        this%setXc => left_flow_channel_interior

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)
      else if (ny == param%ny) then
        this%setXc => left_flow_channel_top_simmetry

        pu = param%hy*(ny + 0.5_DP - param%my2)
        pd = 2.0_DP*(param%ymax - param%hy*(param%my2 - 0.5_DP)) - pu
        this%v(1) = 6.0_DP*param%m*pu*pd
        this%v(2) = 6.0_DP*param%m*pu*pd
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        rr = sqrt(((nx - 1)*param%hx - param%xmax - param%xhs)**2 + &
                  ((ny - 1)*param%hy - (1.0 - param%yhs))**2)
        this%T = param%t0hs*exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP*param%t0hs*exp(-rr/param%r0hs)

      else
        this%setXc => wall_no_exchange

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

        this%T = 0.01_DP
        this%F = 0.0_DP
        this%Z = 0.01_DP
      end if

    end if

    this%setYc => dumb_sub

  end subroutine

  subroutine dumb_sub(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

  end subroutine dumb_sub
  !Left
  subroutine left_flow_channel_interior(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*(point%gammaTx + point%gammaTx) + &
                  param%xy*(point%gammaTy + point%gammaTy)
    point%aF(5) = ap + param%yx*(point%gammaFx + point%gammaFx) + &
                  param%xy*(point%gammaFy + point%gammaFy) + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*(point%gammaZx + point%gammaZx) + &
                  param%xy*(point%gammaZy + point%gammaZy) - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_interior

  subroutine left_flow_channel_inlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = ae - 4.0_DP/3.0_DP*param%yx*point%gammaTx
    point%aF(2) = ae - 4.0_DP/3.0_DP*param%yx*point%gammaFx
    point%aZ(2) = ae - 4.0_DP/3.0_DP*param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*4.0_DP*point%gammaTx + &
                  param%xy*(point%gammaTy + point%gammaTy)
    point%aF(5) = ap + param%yx*4.0_DP*point%gammaFx + &
                  param%xy*(point%gammaFy + point%gammaFy) + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*4.0_DP*point%gammaZx + &
                  param%xy*(point%gammaZy + point%gammaZy) - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 8.0_DP/3.0_DP*point%gammaFx*param%yx + param%hy*param%m*uw
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_inlet

  subroutine left_flow_channel_outlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*point%gammaTx + &
                  param%xy*(point%gammaTy + point%gammaTy)
    point%aF(5) = ap + param%yx*point%gammaFx + &
                  param%xy*(point%gammaFy + point%gammaFy) + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*point%gammaZx + &
                  param%xy*(point%gammaZy + point%gammaZy) - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_outlet

  subroutine left_flow_channel_bot_wall(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*(point%gammaTx + point%gammaTx) + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*(point%gammaFx + point%gammaFx) + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*(point%gammaZx + point%gammaZx) + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_bot_wall

  subroutine left_flow_channel_bot_wall_inlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*4.0_DP*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*4.0_DP*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*4.0_DP*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 8.0_DP/3.0_DP*point%gammaFx*param%yx + param%hy*param%m*uw
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_bot_wall_inlet

  subroutine left_flow_channel_bot_wall_outlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_bot_wall_outlet

  subroutine left_flow_channel_top_simmetry(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = ap + param%yx*(point%gammaTx + point%gammaTx) + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*(point%gammaFx + point%gammaFx) + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*(point%gammaZx + point%gammaZx) + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_top_simmetry

  subroutine left_flow_channel_top_simmetry_inlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = ae - 4.0_DP/3.0_DP*param%yx*point%gammaTx
    point%aF(2) = ae - 4.0_DP/3.0_DP*param%yx*point%gammaFx
    point%aZ(2) = ae - 4.0_DP/3.0_DP*param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = ap + param%yx*4.0_DP*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*4.0_DP*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*4.0_DP*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 8.0_DP/3.0_DP*point%gammaFx*param%yx + param%hy*param%m*uw
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_top_simmetry_inlet

  subroutine left_flow_channel_top_simmetry_outlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = -param%m*param%hy*upwind(uw, "+")
    ae = -param%m*param%hy*upwind(ue, "-")
    as = -param%m*param%hx*upwind(vs, "+")
    an = -param%m*param%hx*upwind(vn, "-")

    ap = param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = ap + param%yx*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine left_flow_channel_top_simmetry_outlet

  ! RIght channel
  subroutine right_flow_channel_interior(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*(point%gammaTx + point%gammaTx) + &
                  param%xy*(point%gammaTy + point%gammaTy)
    point%aF(5) = ap + param%yx*(point%gammaFx + point%gammaFx) + &
                  param%xy*(point%gammaFy + point%gammaFy) + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*(point%gammaZx + point%gammaZx) + &
                  param%xy*(point%gammaZy + point%gammaZy) - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_interior

  subroutine right_flow_channel_inlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - 4.0_DP/3.0_DP*param%yx*point%gammaTx
    point%aF(1) = aw - 4.0_DP/3.0_DP*param%yx*point%gammaFx
    point%aZ(1) = aw - 4.0_DP/3.0_DP*param%yx*point%gammaZx

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*4.0_DP*point%gammaTx + &
                  param%xy*(point%gammaTy + point%gammaTy)
    point%aF(5) = ap + param%yx*4.0_DP*point%gammaFx + &
                  param%xy*(point%gammaFy + point%gammaFy) + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*4.0_DP*point%gammaZx + &
                  param%xy*(point%gammaZy + point%gammaZy) - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 8.0_DP/3.0_DP*point%gammaFx*param%yx + -param%hy*param%m*ue
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_inlet

  subroutine right_flow_channel_outlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = aw - param%yx*point%gammaTx
    point%aF(2) = aw - param%yx*point%gammaFx
    point%aZ(2) = aw - param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*point%gammaTx + &
                  param%xy*(point%gammaTy + point%gammaTy)
    point%aF(5) = ap + param%yx*point%gammaFx + &
                  param%xy*(point%gammaFy + point%gammaFy) + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*point%gammaZx + &
                  param%xy*(point%gammaZy + point%gammaZy) - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_outlet

  subroutine right_flow_channel_bot_simmetry(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*(point%gammaTx + point%gammaTx) + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*(point%gammaFx + point%gammaFx) + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*(point%gammaZx + point%gammaZx) + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_bot_simmetry

  subroutine right_flow_channel_bot_simmetry_inlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*4.0_DP*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*4.0_DP*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*4.0_DP*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 8.0_DP/3.0_DP*point%gammaFx*param%yx - param%hy*param%m*ue
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_bot_simmetry_inlet

  subroutine right_flow_channel_bot_simmetry_outlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = an - param%xy*point%gammaTy
    point%aF(4) = an - param%xy*point%gammaFy
    point%aZ(4) = an - param%xy*point%gammaZy

    point%aT(5) = ap + param%yx*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_bot_simmetry_outlet

  subroutine right_flow_channel_top_wall(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - param%yx*point%gammaTx
    point%aF(1) = aw - param%yx*point%gammaFx
    point%aZ(1) = aw - param%yx*point%gammaZx

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = ap + param%yx*(point%gammaTx + point%gammaTx) + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*(point%gammaFx + point%gammaFx) + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*(point%gammaZx + point%gammaZx) + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_top_wall

  subroutine right_flow_channel_top_wall_inlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         -param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = aw - 4.0_DP/3.0_DP*param%yx*point%gammaTx
    point%aF(1) = aw - 4.0_DP/3.0_DP*param%yx*point%gammaFx
    point%aZ(1) = aw - 4.0_DP/3.0_DP*param%yx*point%gammaZx

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = ap + param%yx*4.0_DP*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*4.0_DP*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*4.0_DP*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 8.0_DP/3.0_DP*point%gammaFx*param%yx - param%hy*param%m*ue
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_top_wall_inlet

  subroutine right_flow_channel_top_wall_outlet(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)

    aw = param%m*param%hy*upwind(uw, "+")
    ae = param%m*param%hy*upwind(ue, "-")
    as = param%m*param%hx*upwind(vs, "+")
    an = param%m*param%hx*upwind(vn, "-")

    ap = -param%m*param%hy*(upwind(uw, "-") + upwind(ue, "+")) + &
         param%m*param%hx*(upwind(vs, "-") + upwind(vn, "+"))

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = ae - param%yx*point%gammaTx
    point%aF(2) = ae - param%yx*point%gammaFx
    point%aZ(2) = ae - param%yx*point%gammaZx

    point%aT(3) = as - param%xy*point%gammaTy
    point%aF(3) = as - param%xy*point%gammaFy
    point%aZ(3) = as - param%xy*point%gammaZy

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = ap + param%yx*point%gammaTx + &
                  param%xy*point%gammaTy
    point%aF(5) = ap + param%yx*point%gammaFx + &
                  param%xy*point%gammaFy + &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%Z*param%hx*param%hy
    point%aZ(5) = ap + param%yx*point%gammaZx + &
                  param%xy*point%gammaZy - &
                  reactionrateFtoZ(param%beta, param%gamma, point%T)*point%F*param%hx*param%hy &
                  + param%hx*param%hy

    point%sT = param%q*point%Z*param%hx*param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine right_flow_channel_top_wall_outlet

  ! Wall
  subroutine wall_no_exchange(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = 1.0_DP
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP

  end subroutine wall_no_exchange

  subroutine wall_exchange(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = -param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = -param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = -param%a2*param%xy
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = -param%a2*param%xy
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = 2.0_DP*param%yx + 2.0_DP*param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange

  subroutine wall_exchange_left_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = -param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = -param%a2*param%xy
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = -param%a2*param%xy
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = param%yx + 2.0_DP*param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_left_boundary

  subroutine wall_exchange_right_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = -param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = -param%a2*param%xy
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = -param%a2*param%xy
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = param%yx + 2.0_DP*param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_right_boundary

  subroutine wall_exchange_top_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = -param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = -param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = -param%a2*param%xy
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = 2.0_DP*param%yx + param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_top_boundary

  subroutine wall_exchange_top_left_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = -param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = -param%a2*param%xy
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = param%yx + param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_top_left_boundary

  subroutine wall_exchange_top_right_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = -param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = -param%a2*param%xy
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = 0.0_DP
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = param%yx + param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_top_right_boundary

  subroutine wall_exchange_bot_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = -param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = -param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = -param%a2*param%xy
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = 2.0_DP*param%yx + param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_bot_boundary

  subroutine wall_exchange_bot_left_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = -param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = -param%a2*param%xy
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = param%yx + param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_bot_left_boundary

  subroutine wall_exchange_bot_right_boundary(point, param)
    class(point_t), intent(inout) :: point
    type(param_t), intent(in) :: param

    point%aT(1) = -param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP

    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP

    point%aT(3) = 0.0_DP
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP

    point%aT(4) = -param%a2*param%xy
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP

    point%aT(5) = param%yx + param%a2*param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP

    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine wall_exchange_bot_right_boundary

  real(DP) function reactionrateFtoZ(beta, gamma, T)
    real(DP) ::beta, gamma, T

    reactionrateFtoZ = beta**2*dexp(beta*(T - 1.0_DP)/ &
                                    (1.0_DP + gamma*(T - 1.0_DP)))
  end function reactionrateFtoZ

  ! Upwind function
  function upwind(v, sign)
    real(DP), intent(in) :: v
    character(len=*), intent(in) :: sign
    real(DP) :: upwind

    if (sign == "+") then
      upwind = (v + abs(v))/2.0_DP
    else if (sign == "-") then
      upwind = (v - abs(v))/2.0_DP
    else
      print *, "How u made this error?"
    end if
  end function upwind

  logical function wall_bot_bd(f)
    procedure(co), pointer :: f
    if (associated(f, wall_exchange_bot_boundary) .or. &
        associated(f, wall_exchange_bot_left_boundary) .or. &
        associated(f, wall_exchange_bot_right_boundary)) then
      wall_bot_bd = .True.
    else
      wall_bot_bd = .False.
    end if
  end function wall_bot_bd

  logical function wall_top_bd(f)
    procedure(co), pointer :: f
    if (associated(f, wall_exchange_top_boundary) .or. &
        associated(f, wall_exchange_top_left_boundary) .or. &
        associated(f, wall_exchange_top_right_boundary)) then
      wall_top_bd = .True.
    else
      wall_top_bd = .False.
    end if
  end function wall_top_bd
end module point2_m
