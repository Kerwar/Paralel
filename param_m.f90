module param_m

  use type_m, only: DP

  implicit none

  private

  type, public :: param_t
    real(DP) :: m, q
    real(DP) :: LeF, LeZ
    real(DP) :: gamma, beta, bK
    real(DP) :: alpha
    real(DP) :: xmax, ymax, a
    real(DP) :: hx, hy, ht

    real(DP) :: palpha, valpha
    real(DP) :: tol
    integer  :: nx, ny, nt
    integer  :: show

    real(DP) :: a2, a4
    ! xy = hx/hy; yx = hy/hx, txy = hx*hy/ht
    real(DP) :: xy, yx, txy
    ! Indeces of exchange in x axis
    integer  :: mx1, mx2
    ! Indeces of wall in y axis
    integer  :: my1, my2

    ! Hot Spot parameters
    real(DP) :: t0hs, xhs, yhs, r0hs

  contains

    procedure :: read_param
  end type param_t

contains

  subroutine read_param(this)
    class(param_t), intent(inout) :: this

    open (101, file="datainput")
    read (101, *)
    read (101, *) this%nx, this%ny, this%nt
    read (101, *)
    read (101, *) this%xmax, this%ymax, this%a
    read (101, *)
    read (101, *) this%m, this%q
    read (101, *)
    read (101, *) this%beta, this%gamma, this%bK
    read (101, *)
    read (101, *) this%LeF, this%LeZ
    read (101, *)
    read (101, *) this%tol, this%show
    read (101, *)
    read (101, *) this%palpha, this%valpha
    read (101, *)
    read (101, *) this%t0hs, this%xhs, this%yhs, this%r0hs
    close (101)

    this%hx = 2.0_DP*this%xmax/(this%nx)
    this%hy = this%ymax/(this%ny)

    this%ht = 10d-6!this%hx * this%hy * 0.00010_DP

    this%a2 = (1.0_DP/this%a)**2
    this%a4 = (1.0_DP/this%a)**4

    this%xy = this%hx/this%hy
    this%yx = this%hy/this%hx
    this%txy = this%hx*this%hy/this%ht

    this%alpha = 0.2_DP
    this%mx1 = int(float(this%nx)/3.0_DP)
    this%mx2 = int(float(this%nx)*2.0_DP/3.0_DP)

    this%my1 = int(float(this%ny)*3.0_DP/7.0_DP)
    this%my2 = this%ny - this%my1!int(float(this%ny) * 4.0_DP/7.0_DP)
  end subroutine read_param

end module param_m
