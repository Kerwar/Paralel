module point_m
  
  use type_m , only: DP
  use param_m, only: param_t
  use paralel_m, only: parll_t
  use interfaces_m, only: coef 
  
  implicit none  
  
  type, public :: point_t
  real(DP) :: T, F, Z
  ! u_w, u_e, u_s, u_n
  real(DP) :: v(4)
  ! W, E, S, N, P
  real(DP) :: aT(5), aF(5), aZ(5), sT, sF, sZ
  
  real(DP) :: gammaTw, gammaTe, gammaTs, gammaTn
  real(DP) :: gammaFw, gammaFe, gammaFs, gammaFn
  real(DP) :: gammaZw, gammaZe, gammaZs, gammaZn
  procedure(coef), nopass, pointer :: cT, cF, cZ
  procedure(co), pointer :: setXc, setYc
  
  contains
  
  procedure :: set_coef
  !procedure :: gamma => init_gamma
  end type point_t
  
  abstract interface
  subroutine co(this, param)
    import :: DP, param_t, point_t, parll_t
    class(point_t), intent(inout) :: this
    type(param_t) , intent(in   ) :: param
  end subroutine
  end interface
  
  contains
  
  subroutine init_gamma(this, param)
    type(point_t), intent(inout) :: this
    type(param_t) , intent(in   ) :: param
    
    this%gammaTw = 1.0_DP
    this%gammaTe = 1.0_DP
    
    this%gammaTs = param%a2
    this%gammaTn = param%a2
    
    this%gammaFw = 1.0_DP/param%LeF
    this%gammaFe = 1.0_DP/param%LeF
    
    this%gammaFs = param%a2/param%LeF
    this%gammaFn = param%a2/param%LeF
    
    this%gammaZw = 1.0_DP/param%LeZ
    this%gammaZe = 1.0_DP/param%LeZ
    
    this%gammaZs = param%a2/param%LeZ
    this%gammaZn = param%a2/param%LeZ
    
  end subroutine init_gamma
  

  ! --------------------------------------------
  !                    ----> 
  ! -------------------------------------------- my2 does not belong to the wall
  !                --------------
  ! -------------------------------------------- my1 does not belongs to the wall
  !                    <----
  ! --------------------------------------------


  subroutine set_coef(this, param, nx, ny)
    class(point_t) , intent(inout) :: this
    type(param_t)  , intent(in   ) :: param
    integer        , intent(in   ) :: nx, ny
    real(DP)       :: rr, pu, pd

    
    if (ny <=  param%my1) then

      pu = param%hy * (ny - 0.5_DP + param%my1)

      pd = 2.0_DP * param%hy * param%my1 - pu 
      if (nx == 1) then
       ! Outlet <--- channels 
        this%setXc => coefx_4

        this%v(1) = - 6.0_DP * param%m * pu * pd 
        this%v(2) = - 6.0_DP * param%m * pu * pd 

        this%T = 0.0_DP 
        this%F = 1.0_DP
        this%Z = 0.0_DP
      else if (nx == param%nx) then
       ! Inlet <--- channels 
        this%setXc => coefx_5

        this%v(1) = - 6.0_DP * param%m * pu * pd 
        this%v(2) = - 6.0_DP * param%m * pu * pd 

        this%T = 0.0_DP
        this%F = 1.0_DP
        this%Z = 0.0_DP
      else
       ! Interior of channel
        this%setXc => coefx_1
        
        this%v(1) = - 6.0_DP * param%m * pu * pd 
        this%v(2) = - 6.0_DP * param%m * pu * pd 

        rr = sqrt(((nx - 1 ) * param%hx - param%xmax + param%xhs)**2 +&
        ((ny - 1 ) * param%hy - param%yhs )**2)
        
        this%T = param%t0hs * exp(-rr/param%r0hs)
        this%F = 1.0_DP
        this%Z = 0.1_DP * param%t0hs * exp(-rr/param%r0hs)
      end if
    else if (ny > param%my1 .and. ny < param%my2) then
      if (nx == param%mx1) then
       ! Wall exchange left bd
        this%setXc => coefx_8

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP

        this%T = 1.0_DP
        this%F = 0.0_DP
        this%Z = 0.0_DP
      else if (nx == param%mx2) then
       ! Wall exchange right bd
        this%setXc => coefx_9

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP

        this%T = 1.0_DP
        this%F = 0.0_DP
        this%Z = 0.0_DP
      else if (nx < param%mx1 .or. nx > param%mx2) then
       ! Wall no exchange
        this%setXc => coefx_6
        
        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP

        this%T = 0.0_DP
        this%F = 0.0_DP
        this%Z = 0.0_DP
      else 
       ! Wall exchange
        this%setXc => coefx_7

        this%v(1) = 0.0_DP
        this%v(2) = 0.0_DP
        
        this%T = 1.0_DP
        this%F = 0.0_DP
        this%Z = 0.0_DP
      end if
    else if (ny >= param%my2 ) then
      pu = param%hy * (ny + 0.5_DP - param%my2)

      pd = 2.0_DP * (param%ymax - param%hy * param%my2) - pu
      if (nx == 1) then
       ! Inlet ---> channels
        this%setXc => coefx_2
        this%v(1) = 6.0_DP * param%m * pu * pd 
        this%v(2) = 6.0_DP * param%m * pu * pd 

        this%T = 0.0_DP
        this%F = 1.0_DP
        this%Z = 0.0_DP
      else if (nx == param%nx) then
       ! Outlet ---> channels 
        this%setXc => coefx_3

        this%v(1) = 6.0_DP * param%m * pu * pd 
        this%v(2) = 6.0_DP * param%m * pu * pd 

        this%T = 0.0_DP
        this%F = 1.0_DP
        this%Z = 0.0_DP
      else
       ! Interior of channel
        this%setXc => coefx_1
        
        this%v(1) = 6.0_DP * param%m * pu * pd 
        this%v(2) = 6.0_DP * param%m * pu * pd 

        rr = sqrt(((nx - 1 ) * param%hx - param%xmax - param%xhs)**2 +&
        ((ny - 1 ) * param%hy  -  (1.0 - param%yhs) )**2)
        this%T = param%t0hs * exp(-rr/param%r0hs)
        if (this%T == 0.0_DP) print *, nx, ny
        this%F = 1.0_DP
        this%Z = 0.1_DP * param%t0hs * exp(-rr/param%r0hs)
      end if
      
    end if
    
    if (nx < param%mx1) then
     ! Interior of bot channel
      if (ny < param%my1) then
        this%setYc => coefy_1

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Boundary with wall no exchange
      else if (ny == param%my1 ) then
        this%setYc => coefy_2

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Interior of top channel
      else if (ny > param%my2) then
        this%setYc => coefy_1

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Boundary with wall no exchange
      else if (ny == param%my2) then
        this%setYc => coefy_3

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Wall
      else
        this%setYc => coefy_5

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
      end if
    elseif (nx > param%mx2) then
     ! Interior of bot channel
      if (ny < param%my1) then
        this%setYc => coefy_1

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP

     ! Boundary with wall no exchange
      else if (ny == param%my1) then
        this%setYc => coefy_2

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Interior of top channel
      else if (ny > param%my2) then
        this%setYc => coefy_1

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Boundary with wall no exchange
      else if (ny == param%my2) then
        this%setYc => coefy_3
        
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Wall
      else
        this%setYc => coefy_5

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
      end if
    else 
     ! Bot channel
      if (ny < param%my1) then
        this%setYc => coefy_1
        
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Bot channel top boundary
      else if (ny == param%my1) then
        this%setYc => coefy_6
        
        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Top channel
      else if (ny > param%my2) then
        this%setYc => coefy_1

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Top channel bot boundary
      else if (ny == param%my2) then
        this%setYc => coefy_7

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
     ! Wall
      else
        this%setYc => coefy_4

        this%v(3) = 0.0_DP
        this%v(4) = 0.0_DP
      end if
    end if

  end subroutine
      
 ! Interior of channel
  subroutine coefx_1(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, ap
    
    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)
    
    aw = - param%m * param%hy * upwind(uw, "+")
    ae = - param%m * param%hy * upwind(ue, "-") 
    
    ap = param%m * param%hy * (upwind(uw, "-") + upwind(ue, "+")) + param%txy
    
    point%aT(1) = aw - param%yx * point%gammaTw
    point%aF(1) = aw - param%yx * point%gammaFw
    point%aZ(1) = aw - param%yx * point%gammaZw
    
    point%aT(2) = ae - param%yx * point%gammaTe
    point%aF(2) = ae - param%yx * point%gammaFe
    point%aZ(2) = ae - param%yx * point%gammaZe
    
    point%aT(5) = ap + param%yx * (point%gammaTw + point%gammaTe)
    point%aF(5) = ap + param%yx * (point%gammaFw + point%gammaFe) + &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%Z * param%hx * param%hy
    point%aZ(5) = ap + param%yx * (point%gammaZw + point%gammaZe) - &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%F * param%hx * param%hy &
    + param%hx * param%hy
    
    point%sT = param%q * point%Z * param%hx * param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefx_1
      
 ! Inlet ---> channels 
  subroutine coefx_2(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, ap
    
    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)
    
    aw = 0.0_DP
    ae = - param%m * param%hy * upwind(ue, "-") 
    
    
    ap = param%m * param%hy * upwind(ue, "+") + param%txy
    
    point%aT(1) = aw
    point%aF(1) = aw
    point%aZ(1) = aw
    
    point%aT(2) = ae - param%yx * &
    (point%gammaTe + 1.0_DP/3.0_DP * point%gammaTw)
    point%aF(2) = ae - param%yx * &
    (point%gammaFe + 1.0_DP/3.0_DP * point%gammaFw)
    point%aZ(2) = ae - param%yx * &
    (point%gammaZe + 1.0_DP/3.0_DP * point%gammaZw)
    
    point%aT(5) = ap + param%yx * &
    (3.0_DP * point%gammaTw + point%gammaTe)
    point%aF(5) = ap + param%yx * &
    (3.0_DP * point%gammaFw + point%gammaFe)+ &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%Z * param%hx * param%hy
    point%aZ(5) = ap + param%yx * &
    (3.0_DP * point%gammaZw + point%gammaZe) - &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%F * param%hx * param%hy &
    + param%hx * param%hy
    
    point%sT = param%q * point%Z * param%hx * param%hy
    point%sF = 8.0_DP/3.0_DP * point%gammaFw * param%yx - param%hy * param%m * point%v(1) 
    point%sZ = 0.0_DP
  end subroutine coefx_2
      
 ! Outlet ---> channels 
  subroutine coefx_3(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, ap
    
    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)
    
    aw = - param%m * param%hy * upwind(uw, "+") 
    ae = 0.0_DP
    
    ap = param%m * param%hy * (upwind(uw, "-") + ue) + param%txy
    
    point%aT(1) = aw - param%yx * point%gammaTw
    point%aF(1) = aw - param%yx * point%gammaFw
    point%aZ(1) = aw - param%yx * point%gammaZw
    
    point%aT(2) = ae
    point%aF(2) = ae
    point%aZ(2) = ae
    
    point%aT(5) = ap + param%yx * (point%gammaTw)
    point%aF(5) = ap + param%yx * (point%gammaFw)+ &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%Z * param%hx * param%hy
    point%aZ(5) = ap + param%yx * (point%gammaZw) - &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%F * param%hx * param%hy &
    + param%hx * param%hy
    
    point%sT = param%q * point%Z * param%hx * param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefx_3
      
 ! Outlet <--- channels 
  subroutine coefx_4(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, ap
    
    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)
    
    aw = 0.0_DP
    ae = - param%m * param%hy * upwind(uw, "-") 
    
    ap = param%m * param%hy * (upwind(ue, "+") + uw) + param%txy
    
    point%aT(1) = aw
    point%aF(1) = aw
    point%aZ(1) = aw

    point%aT(2) = ae - param%yx * point%gammaTe
    point%aF(2) = ae - param%yx * point%gammaFe
    point%aZ(2) = ae - param%yx * point%gammaZe
    
    
    point%aT(5) = ap + param%yx * (point%gammaTe)
    point%aF(5) = ap + param%yx * (point%gammaFe) + &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%Z * param%hx * param%hy
    point%aZ(5) = ap + param%yx * (point%gammaZe) - &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%F * param%hx * param%hy &
    + param%hx * param%hy
    
    point%sT = param%q * point%Z * param%hx * param%hy
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefx_4
      
 ! Inlet <--- channels 
  subroutine coefx_5(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: uw, ue, vs, vn
    real(DP) :: aw, ae, ap
    
    uw = point%v(1)
    ue = point%v(2)
    vs = point%v(3)
    vn = point%v(4)
    
    aw = - param%m * param%hy * upwind(uw, "+") 
    ae = 0.0_DP
    
    ap = param%m * param%hy * upwind(uw, "-") + param%txy
    
    point%aT(1) = aw - param%yx * &
    (point%gammaTw + 1.0_DP/3.0_DP * point%gammaTe)
    point%aF(1) = aw - param%yx * &
    (point%gammaFw + 1.0_DP/3.0_DP * point%gammaFe) 
    point%aZ(1) = aw - param%yx * &
    (point%gammaZw + 1.0_DP/3.0_DP * point%gammaZe)
    
    point%aT(2) = ae
    point%aF(2) = ae
    point%aZ(2) = ae
    
    point%aT(5) = ap + param%yx * &
    (3.0_DP * point%gammaTe + point%gammaTw)
    point%aF(5) = ap + param%yx * &
    (3.0_DP * point%gammaFe + point%gammaFw)+ &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%Z * param%hx * param%hy
    point%aZ(5) = ap + param%yx * &
    (3.0_DP * point%gammaZe + point%gammaZw) - &
    param%beta**2 * dexp(param%beta * (point%T - 1.0_DP) /          &
    (1.0_DP + param%gamma * (point%T - 1.0_DP))) * point%F * param%hx * param%hy &
    + param%hx * param%hy
    
    point%sT = param%q * point%Z * param%hx * param%hy
    point%sF = 8.0_DP/3.0_DP * point%gammaFe * param%yx - param%hy * param%m * point%v(1)
    point%sZ = 0.0_DP
  end subroutine coefx_5
  
 ! Wall no exchange
  subroutine coefx_6(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP
    
    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP
    
    point%aT(5) = 1.0_DP
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP
    
    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefx_6
  
 ! Wall exchange
  subroutine coefx_7(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    
    point%aT(1) = - param%alpha * param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP
    
    point%aT(2) = - param%alpha * param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP
    
    point%aT(5) = 2.0_DP * param%alpha * param%yx + param%txy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP
    
    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefx_7
  
 ! Wall exchange left bd
  subroutine coefx_8(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    
    point%aT(1) = 0.0_DP
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP
    
    point%aT(2) = - (4.0_DP/3.0_DP) * param%alpha * param%yx
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP
    
    point%aT(5) = 4.0_DP * param%alpha * param%yx + param%txy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP
    
    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefx_8
  
 ! Wall exchange right bd
  subroutine coefx_9(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    
    point%aT(1) = - (4.0_DP/3.0_DP) * param%alpha * param%yx
    point%aF(1) = 0.0_DP
    point%aZ(1) = 0.0_DP
    
    point%aT(2) = 0.0_DP
    point%aF(2) = 0.0_DP
    point%aZ(2) = 0.0_DP
    
    point%aT(5) = 4.0_DP * param%alpha * param%yx + param%txy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP
    
    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefx_9
      
 ! Interior of channel
  subroutine coefy_1(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: vs, vn
    real(DP) :: as, an, ap
    
    vs = point%v(3)
    vn = point%v(4)
    
    as = - param%m * param%hx * upwind(vs, "+")
    an = - param%m * param%hx * upwind(vn, "-") 
    
    ap = param%m * param%hx * (upwind(vs, "-") + upwind(vn, "+")) 
    
    point%aT(3) = as - param%xy * point%gammaTs
    point%aF(3) = as - param%xy * point%gammaFs
    point%aZ(3) = as - param%xy * point%gammaZs
    
    point%aT(4) = an - param%xy * point%gammaTn
    point%aF(4) = an - param%xy * point%gammaFn
    point%aZ(4) = an - param%xy * point%gammaZn
    
    point%aT(5) = point%aT(5) + ap + &
    param%xy * (point%gammaTs + point%gammaTn)
    point%aF(5) = point%aF(5) + ap + &
    param%xy * (point%gammaFs + point%gammaFn)
    point%aZ(5) = point%aZ(5) + ap + &
    param%xy * (point%gammaZs + point%gammaZn)
    
    point%sT = point%sT + 0.0_DP
    point%sF = point%sF + 0.0_DP
    point%sZ = point%sZ + 0.0_DP
  end subroutine coefy_1
  
 ! Top bd no exchange
  subroutine coefy_2(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: vs, vn
    real(DP) :: as, an, ap
    
    vs = point%v(3)
    vn = point%v(4)
    
    as = - param%m * param%hx * upwind(vs, "+")
    an = 0.0_DP
    
    ap = param%m * param%hx * upwind(vs, "-") 
    
    point%aT(3) = as - param%xy * &
    (point%gammaTs + 1.0_DP/3.0_DP * point%gammaTn)
    point%aF(3) = as - param%xy * &
    (point%gammaFs                                )
    point%aZ(3) = as - param%xy * & 
    (point%gammaZs                                )
    
    point%aT(4) = an 
    point%aF(4) = an 
    point%aZ(4) = an 
    
    point%aT(5) = point%aT(5) + ap + &
    param%xy * (point%gammaTs + 3.0_DP * point%gammaTn)
    point%aF(5) = point%aF(5) + ap + &
    param%xy * (point%gammaFs                         )
    point%aZ(5) = point%aZ(5) + ap + &
    param%xy * (point%gammaZs                         )
    
    point%sT = point%sT + 0.0_DP
    point%sF = point%sF + 0.0_DP
    point%sZ = point%sZ + 0.0_DP
  end subroutine coefy_2
  
 ! Bot bd no exchange
  subroutine coefy_3(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: vs, vn
    real(DP) :: as, an, ap
    
    vs = point%v(3)
    vn = point%v(4)
    
    as = 0.0_DP
    an = - param%m * param%hx * upwind(vn, "-")
    
    ap = param%m * param%hx * upwind(vn, "+") 
    
    point%aT(3) = as 
    point%aF(3) = as 
    point%aZ(3) = as 
    
    point%aT(4) = an - param%xy * &
    (1.0_DP/3.0_DP * point%gammaTs + point%gammaTn)
    point%aF(4) = an - param%xy * &
    (                                point%gammaFn)
    point%aZ(4) = an - param%xy * &
    (                                point%gammaZn)
    
    point%aT(5) = point%aT(5) + ap + &
    param%xy * (3.0_DP * point%gammaTs + point%gammaTn)
    point%aF(5) = point%aF(5) + ap + &
    param%xy * (                         point%gammaFn)
    point%aZ(5) = point%aZ(5) + ap + &
    param%xy * (                         point%gammaZn)
    
    point%sT = point%sT + 0.0_DP
    point%sF = point%sF + 0.0_DP
    point%sZ = point%sZ + 0.0_DP
  end subroutine coefy_3
      
 ! Wall exchange
  subroutine coefy_4(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    
    point%aT(3) = - param%alpha * param%a2 * param%xy
    point%aF(3) = 0.0_DP
    point%aZ(3) = 0.0_DP
    
    point%aT(4) = - param%alpha * param%a2 * param%xy
    point%aF(4) = 0.0_DP
    point%aZ(4) = 0.0_DP
    
    point%aT(5) = point%aT(5) + 2.0_DP * param%a2 * param%alpha * param%xy
    point%aF(5) = 1.0_DP
    point%aZ(5) = 1.0_DP
    
    point%sT = 0.0_DP
    point%sF = 0.0_DP
    point%sZ = 0.0_DP
  end subroutine coefy_4
  
 ! Wall no exchange
  subroutine coefy_5(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    
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
  end subroutine coefy_5

 ! Top bd exchange
  subroutine coefy_6(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: vs, vn
    real(DP) :: as, an, ap
    
    vs = point%v(3)
    vn = point%v(4)
    
    as = - param%m * param%hx * upwind(vs, "+")
    an = 0.0_DP
    
    ap = param%m * param%hx * upwind(vs, "-") 
    
    point%aT(3) = as - param%xy * &
    (point%gammaTs + 1.0_DP/3.0_DP * point%gammaTn)
    point%aF(3) = as - param%xy * &
    (point%gammaFs                                )
    point%aZ(3) = as - param%xy * & 
    (point%gammaZs                                )
    
    point%aT(4) = an 
    point%aF(4) = an 
    point%aZ(4) = an 
    
    point%aT(5) = point%aT(5) + ap + &
    param%xy * (point%gammaTs + 3.0_DP * point%gammaTn)
    point%aF(5) = point%aF(5) + ap + &
    param%xy * (point%gammaFs                         )
    point%aZ(5) = point%aZ(5) + ap + &
    param%xy * (point%gammaZs                         )
    
    point%sT = point%sT + 0.0_DP
    point%sF = point%sF + 0.0_DP
    point%sZ = point%sZ + 0.0_DP
  end subroutine coefy_6

 ! Bot bd exchange
  subroutine coefy_7(point, param)
    class(point_t), intent(inout) :: point
    type(param_t) , intent(in   ) :: param
    
    real(DP) :: vs, vn
    real(DP) :: as, an, ap
    
    vs = point%v(3)
    vn = point%v(4)
    
    as = 0.0_DP
    an = - param%m * param%hx * upwind(vn, "-")
    
    ap = param%m * param%hx * upwind(vn, "+") 
    
    point%aT(3) = as 
    point%aF(3) = as 
    point%aZ(3) = as 
    
    point%aT(4) = an - param%xy * &
    (1.0_DP/3.0_DP * point%gammaTs + point%gammaTn)
    point%aF(4) = an - param%xy * &
    (                                point%gammaFn)
    point%aZ(4) = an - param%xy * &
    (                                 point%gammaZn)
    
    point%aT(5) = point%aT(5) + ap + &
    param%xy * (3.0_DP * point%gammaTs + point%gammaTn)
    point%aF(5) = point%aF(5) + ap + &
    param%xy * (                         point%gammaFn)
    point%aZ(5) = point%aZ(5) + ap + &
    param%xy * (                         point%gammaZn)
    
    point%sT = point%sT + 0.0_DP
    point%sF = point%sF + 0.0_DP  
    point%sZ = point%sZ + 0.0_DP
  end subroutine coefy_7
 ! Upwind function
  function upwind(v, sign)
    real(DP), intent(in   ) :: v
    character(len = *), intent(in   ) :: sign
    real(DP) :: upwind
    
    if (sign == "+") then
      upwind = (v + abs(v)) / 2.0_DP
    else if (sign == "-") then
      upwind = (v - abs(v)) / 2.0_DP
    else
      print *, "How u made this error?"
    end if
  end function upwind
end module point_m