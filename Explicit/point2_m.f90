module point_m

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

  contains

    procedure :: set_coef
    procedure, private :: set_to_left_channel, set_to_right_channel, set_to_wall_channel
    !procedure :: gamma => init_gamma
  end type point_t

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

    if (ny >= param%my2) then
      call this%set_to_right_channel(param, nx, ny)
    else if (ny < param%my1) then
      call this%set_to_left_channel(param, nx, ny)
    else
      call this%set_to_wall_channel(param, nx, ny)
    end if

  end subroutine

  subroutine set_to_right_channel(this, param, nx, ny)
    class(point_t), intent(inout) :: this
    type(param_t), intent(in) :: param
    integer, intent(in) :: nx, ny

    real(DP) :: m, gammaTx, gammaFx, gammaZx, gammaTy, gammaFy, gammaZy, beta, hx, hy
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = this%v(1)
    ue = this%v(2)
    vs = this%v(3)
    vn = this%v(4)

    m = param%m
    gammaTx = this%gammaTx
    gammaFx = this%gammaFx
    gammaZx = this%gammaZx
    gammaTy = this%gammaTy
    gammaFy = this%gammaFy
    gammaZy = this%gammaZy
    beta = param%beta
    hx = param%hx
    hy = param%hy

    if (ny == param%ny) then

    else if (ny == param%my2) then
      if (nx >= param%mx1 .and. nx < param%mx2) then
      else
      end if
    else if (nx == 1) then
      aw = 0.0_DP
      ae = -m*hy*upwind(ue, "-")
      as = -m*hx*upwind(vs, "+")
      an = -m*hx*upwind(vn, "-")

      ap = -(aw + ae + as + an) + param%txy

      this%aT(1) = aw
      this%aF(1) = aw
      this%aZ(1) = aw

      this%aT(2) = ae - param%yx* &
                   (gammaTx + 1.0_DP/3.0_DP*gammaTx)
      this%aF(2) = ae - param%yx* &
                   (gammaTx + 1.0_DP/3.0_DP*gammaTx)
      this%aZ(2) = ae - param%yx* &
                   (gammaTx + 1.0_DP/3.0_DP*gammaTx)

      this%aT(3) = as - param%xy*gammaTy
      this%aF(3) = as - param%xy*gammaFy
      this%aZ(3) = as - param%xy*gammaZy

      this%aT(4) = an - param%xy*gammaTy
      this%aF(4) = an - param%xy*gammaFy
      this%aZ(4) = an - param%xy*gammaZy

      this%aT(5) = ap + param%yx*(3.0_DP*gammaTx + gammaTx) + &
                   param%xy*(gammaTy + gammaTy)
      this%aF(5) = ap + param%yx*(3.0_DP*gammaFx + gammaFx) + &
                   param%xy*(gammaFy + gammaFy) + &
                   beta**2*dexp(beta*(this%T - 1.0_DP)/ &
                                (1.0_DP + param%gamma*(this%T - 1.0_DP)))*this%Z*hx*hy
      this%aZ(5) = ap + param%yx*(3.0_DP*gammaZx + gammaZx) + &
                   param%xy*(gammaZy + gammaZy) - &
                   beta**2*dexp(beta*(this%T - 1.0_DP)/ &
                                (1.0_DP + param%gamma*(this%T - 1.0_DP)))*this%F*hx*hy &
                   + hx*hy

      this%sT = param%q*this%Z*hx*hy
      this%sF = 8.0_DP/3.0_DP*gammaFx*param%yx + param%hy*param%m*uw
      this%sZ = 0.0_DP
    else if (nx == param%nx) then
      aw = -m*hy*upwind(uw, "+")
      ae = -m*hy*ue
      as = -m*hx*upwind(vs, "+")
      an = -m*hx*upwind(vn, "-")

      ap = -(aw + ae + as + an) + param%txy

      this%aT(1) = aw - param%yx*gammaTx
      this%aF(1) = aw - param%yx*gammaFx
      this%aZ(1) = aw - param%yx*gammaZx

      this%aT(2) = 0.0_DP
      this%aF(2) = 0.0_DP
      this%aZ(2) = 0.0_DP

      this%aT(3) = as - param%xy*gammaTy
      this%aF(3) = as - param%xy*gammaFy
      this%aZ(3) = as - param%xy*gammaZy

      this%aT(4) = an - param%xy*gammaTy
      this%aF(4) = an - param%xy*gammaFy
      this%aZ(4) = an - param%xy*gammaZy

      this%aT(5) = ap + param%yx*(gammaTx) + &
                   param%xy*(gammaTy + gammaTy)
      this%aF(5) = ap + param%yx*(gammaFx) + &
                   param%xy*(gammaFy + gammaFy) + &
                   beta**2*dexp(beta*(this%T - 1.0_DP)/ &
                                (1.0_DP + param%gamma*(this%T - 1.0_DP)))*this%Z*hx*hy
      this%aZ(5) = ap + param%yx*(gammaZx) + &
                   param%xy*(gammaZy + gammaZy) - &
                   beta**2*dexp(beta*(this%T - 1.0_DP)/ &
                                (1.0_DP + param%gamma*(this%T - 1.0_DP)))*this%F*hx*hy &
                   + hx*hy

      this%sT = param%q*this%Z*hx*hy
      this%sF = 0.0_DP
      this%sZ = 0.0_DP
    else ! Interior of channel
      aw = -m*hy*upwind(uw, "+")
      ae = -m*hy*upwind(ue, "-")
      as = -m*hx*upwind(vs, "+")
      an = -m*hx*upwind(vn, "-")

      ap = -(aw + ae + as + an) + param%txy

      this%aT(1) = aw - param%yx*gammaTx
      this%aF(1) = aw - param%yx*gammaFx
      this%aZ(1) = aw - param%yx*gammaZx

      this%aT(2) = ae - param%yx*gammaTx
      this%aF(2) = ae - param%yx*gammaFx
      this%aZ(2) = ae - param%yx*gammaZx

      this%aT(3) = as - param%xy*gammaTy
      this%aF(3) = as - param%xy*gammaFy
      this%aZ(3) = as - param%xy*gammaZy

      this%aT(4) = an - param%xy*gammaTy
      this%aF(4) = an - param%xy*gammaFy
      this%aZ(4) = an - param%xy*gammaZy

      this%aT(5) = ap + param%yx*(gammaTx + gammaTx) + &
                   param%xy*(gammaTy + gammaTy)
      this%aF(5) = ap + param%yx*(gammaFx + gammaFx) + &
                   param%xy*(gammaFy + gammaFy) + &
                   beta**2*dexp(beta*(this%T - 1.0_DP)/ &
                                (1.0_DP + param%gamma*(this%T - 1.0_DP)))*this%Z*hx*hy
      this%aZ(5) = ap + param%yx*(gammaZx + gammaZx) + &
                   param%xy*(gammaZy + gammaZy) - &
                   beta**2*dexp(beta*(this%T - 1.0_DP)/ &
                                (1.0_DP + param%gamma*(this%T - 1.0_DP)))*this%F*hx*hy &
                   + hx*hy

      this%sT = param%q*this%Z*hx*hy
      this%sF = 0.0_DP
      this%sZ = 0.0_DP
    end if
  end subroutine set_to_right_channel

  subroutine set_to_left_channel(this, param, nx, ny)
    class(point_t), intent(inout) :: this
    type(param_t), intent(in) :: param
    integer, intent(in) :: nx, ny

    real(DP) :: m, gammaTx, gammaFx, gammaZx, gammaTy, gammaFy, gammaZy, beta, hx, hy
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = this%v(1)
    ue = this%v(2)
    vs = this%v(3)
    vn = this%v(4)

    m = -param%m
    gammaTx = this%gammaTx
    gammaFx = this%gammaFx
    gammaZx = this%gammaZx
    gammaTy = this%gammaTy
    gammaFy = this%gammaFy
    gammaZy = this%gammaZy
    beta = param%beta
    hx = param%hx
    hy = param%hy

    if (ny == 1) then

    else if (ny == param%my1) then
      if (nx >= param%mx1 .and. nx < param%mx2) then
      else
      end if
    else if (nx == param%nx) then
      aw = -m*hy*upwind(uw, "+")
      ae = 0.0_DP
      as = -m*hx*upwind(vs, "+")
      an = -m*hx*upwind(vn, "-")

      ap = -(aw + ae + as + an) + param%txy

      this%aT(1) = aw - param%yx* &
                   (gammaTx + 1.0_DP/3.0_DP*gammaTx)
      this%aF(1) = aw - param%yx* &
                   (gammaTx + 1.0_DP/3.0_DP*gammaTx)
      this%aZ(1) = aw - param%yx* &
                   (gammaTx + 1.0_DP/3.0_DP*gammaTx)

      this%aT(2) = ae
      this%aF(2) = ae
      this%aZ(2) = ae

      this%aT(3) = as - param%xy*gammaTy
      this%aF(3) = as - param%xy*gammaFy
      this%aZ(3) = as - param%xy*gammaZy

      this%aT(4) = an - param%xy*gammaTy
      this%aF(4) = an - param%xy*gammaFy
      this%aZ(4) = an - param%xy*gammaZy

      this%aT(5) = ap + param%yx*(3.0_DP*gammaTx + gammaTx) + &
                   param%xy*(gammaTy + gammaTy)
      this%aF(5) = ap + param%yx*(3.0_DP*gammaFx + gammaFx) + &
                   param%xy*(gammaFy + gammaFy) + &
                   reactionrateFtoZ(beta, param%gamma, this%T)*this%Z*hx*hy
      this%aZ(5) = ap + param%yx*(3.0_DP*gammaZx + gammaZx) + &
                   param%xy*(gammaZy + gammaZy) - &
                   reactionrateFtoZ(beta, param%gamma, this%T)*this%F*hx*hy &
                   + hx*hy

      this%sT = param%q*this%Z*hx*hy
      this%sF = 8.0_DP/3.0_DP*gammaFx*param%yx + param%hy*param%m*uw
      this%sZ = 0.0_DP
    else if (nx == 1) then
      aw = -m*hy*upwind(uw, "+")
      ae = -m*hy*ue
      as = -m*hx*upwind(vs, "+")
      an = -m*hx*upwind(vn, "-")

      ap = -(aw + ae + as + an) + param%txy

      this%aT(1) = aw - param%yx*gammaTx
      this%aF(1) = aw - param%yx*gammaFx
      this%aZ(1) = aw - param%yx*gammaZx

      this%aT(2) = 0.0_DP
      this%aF(2) = 0.0_DP
      this%aZ(2) = 0.0_DP

      this%aT(3) = as - param%xy*gammaTy
      this%aF(3) = as - param%xy*gammaFy
      this%aZ(3) = as - param%xy*gammaZy

      this%aT(4) = an - param%xy*gammaTy
      this%aF(4) = an - param%xy*gammaFy
      this%aZ(4) = an - param%xy*gammaZy

      this%aT(5) = ap + param%yx*(gammaTx) + &
                   param%xy*(gammaTy + gammaTy)
      this%aF(5) = ap + param%yx*(gammaFx) + &
                   param%xy*(gammaFy + gammaFy) + &
                   reactionrateFtoZ(beta, param%gamma, this%T)*this%Z*hx*hy
      this%aZ(5) = ap + param%yx*(gammaZx) + &
                   param%xy*(gammaZy + gammaZy) - &
                   reactionrateFtoZ(beta, param%gamma, this%T)*this%F*hx*hy &
                   + hx*hy

      this%sT = param%q*this%Z*hx*hy
      this%sF = 0.0_DP
      this%sZ = 0.0_DP
    else ! Interior of channel
      aw = -m*hy*upwind(uw, "+")
      ae = -m*hy*upwind(ue, "-")
      as = -m*hx*upwind(vs, "+")
      an = -m*hx*upwind(vn, "-")

      ap = -(aw + ae + as + an) + param%txy

      this%aT(1) = aw - param%yx*gammaTx
      this%aF(1) = aw - param%yx*gammaFx
      this%aZ(1) = aw - param%yx*gammaZx

      this%aT(2) = ae - param%yx*gammaTx
      this%aF(2) = ae - param%yx*gammaFx
      this%aZ(2) = ae - param%yx*gammaZx

      this%aT(3) = as - param%xy*gammaTy
      this%aF(3) = as - param%xy*gammaFy
      this%aZ(3) = as - param%xy*gammaZy

      this%aT(4) = an - param%xy*gammaTy
      this%aF(4) = an - param%xy*gammaFy
      this%aZ(4) = an - param%xy*gammaZy

      this%aT(5) = ap + param%yx*(gammaTx + gammaTx) + &
                   param%xy*(gammaTy + gammaTy)
      this%aF(5) = ap + param%yx*(gammaFx + gammaFx) + &
                   param%xy*(gammaFy + gammaFy) + &
                   reactionrateFtoZ(beta, param%gamma, this%T)*this%Z*hx*hy
      this%aZ(5) = ap + param%yx*(gammaZx + gammaZx) + &
                   param%xy*(gammaZy + gammaZy) - &
                   reactionrateFtoZ(beta, param%gamma, this%T)*this%F*hx*hy &
                   + hx*hy

      this%sT = param%q*this%Z*hx*hy
      this%sF = 0.0_DP
      this%sZ = 0.0_DP
    end if
  end subroutine set_to_left_channel

  subroutine set_to_wall_channel(this, param, nx, ny)
    class(point_t), intent(inout) :: this
    type(param_t), intent(in) :: param
    integer, intent(in) :: nx, ny

    real(DP) :: m, gammaTx, gammaFx, gammaZx, gammaTy, gammaFy, gammaZy, beta, hx, hy
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, as, an, ap

    uw = this%v(1)
    ue = this%v(2)
    vs = this%v(3)
    vn = this%v(4)

    m = param%m
    gammaTx = this%gammaTx
    gammaFx = this%gammaFx
    gammaZx = this%gammaZx
    gammaTy = this%gammaTy
    gammaFy = this%gammaFy
    gammaZy = this%gammaZy
    beta = param%beta
    hx = param%hx
    hy = param%hy

  end subroutine set_to_wall_channel

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

  real(DP) function diffusion(gamma, boundarytype)
    real(DP) :: gamma
    character(*), optional :: boundarytype

    if (present(boundarytype)) then
    else
      diffusion = gamma
    end if
  end function
  real(DP) function reactionrateFtoZ(beta, gamma, T)
    real(DP) ::beta, gamma, T

    reactionrateFtoZ = beta**2*dexp(beta*(T - 1.0_DP)/ &
                                    (1.0_DP + gamma*(T - 1.0_DP)))
  end function reactionrateFtoZ
end module point_m
