module coefun_m

  use type_m, only: DP
  use param_m, only: param_t
  use point_m

  !include 'mpif.h'
  implicit none

contains
  ! Channel 1
  ! 1  -> top-left
  ! 2  -> top-right
  ! 3  -> bot-left
  ! 4  -> bot-right
  ! 5  -> top
  ! 6  -> bot-noex
  ! 7  -> bot-noex-ex
  ! 8  -> bot-ex
  ! 9  -> bot-ex-noex
  ! 10 -> left
  ! 11 -> right
  ! Channel 2
  ! 12 -> top-left
  ! 13 -> top-right
  ! 14 -> bot-left
  ! 15 -> bot-right
  ! 16 -> top-noex
  ! 17 -> top-noex-ex
  ! 18 -> top-ex
  ! 19 -> top-ex-noex
  ! 20 -> bot
  ! 21 -> left
  ! 22 -> right

  ! 23 -> wall left
  ! 24 -> wall right
  ! 25 -> main
  ! --------------------------------------------
  !                    ---->
  ! --------------------------------------------
  !                --------------
  ! --------------------------------------------
  !                    <----
  ! --------------------------------------------

  !subroutine solve(field, PL)
  !  class(point_t), intent(inout) :: field(:,:)
  !  type(parll_t) , intent(in   ) :: PL
  !
  !
  !end subroutine solve

  subroutine ADI(w_d, e_d, s_d, n_d, o_d, ind, sol, error, swps, rlx)
    real(DP), intent(in) :: w_d(:, :), e_d(:, :), o_d(:, :), &
                            s_d(:, :), n_d(:, :), ind(:, :), rlx
    real(DP), intent(inout) :: sol(:, :), error
    integer, intent(in) :: swps
    real(DP), dimension(size(o_d, 1)*size(o_d, 2)) :: a, b, c, d, e, S, x

    integer :: nx, ny, i, j, k, id

    nx = size(o_d, 1)
    ny = size(o_d, 2)

    do k = 1, swps
      do j = 1, ny
        do i = 1, nx

          id = (j - 1)*nx + i
          a(id) = s_d(i, j)
          b(id) = w_d(i, j)
          c(id) = o_d(i, j)
          d(id) = e_d(i, j)
          e(id) = n_d(i, j)

          S(id) = ind(i, j) - o_d(i, j)*sol(i, j)
          ! if (k == 1 ) then
          !   error_v(id) = ind(i,j)-o_d(i,j)*sol(i,j)
          !   if(i > 1)  error_V(id) = error_v(id) - w_d(i,j)*sol(i-1,j)
          !   if(i < nx) error_V(id) = error_v(id) - e_d(i,j)*sol(i+1,j)
          !   if(j > 1)  error_V(id) = error_v(id) - s_d(i,j)*sol(i,j-1)
          !   if(j < ny) error_V(id) = error_v(id) - n_d(i,j)*sol(i,j+1)
          ! end if
          if (i > 1) S(id) = S(id) - w_d(i, j)*sol(i - 1, j)
          if (i < nx) S(id) = S(id) - e_d(i, j)*sol(i + 1, j)
          if (j > 1) S(id) = S(id) - s_d(i, j)*sol(i, j - 1)
          if (j < ny) S(id) = S(id) - n_d(i, j)*sol(i, j + 1)

          !print "(I4, 6ES25.14)", id, a(id), b(id), c(id), d(id), e(id), S(id)
        end do
      end do
      !stop
      if (k == 1) error = norm2(S)!norm2(error_v)

      x = 0.0_DP
      ! do i = 1, nx
      !   do j = 1, ny
      !     x((j-1)*nx + i) = sol(i,j)
      !   end do
      ! end do

      do j = 1, ny

        id = (j - 1)*nx + 1

        if (j > 1) S(id:id + nx - 1) = S(id:id + nx - 1) - &
                                       a(id:id + nx - 1)*x(id - nx:id - 1)
        if (j < ny) S(id:id + nx - 1) = S(id:id + nx - 1) - &
                                        e(id:id + nx - 1)*x(id + nx:id + 2*nx - 1)

        !print *, o_d(:,1)

        call trisolver(b(id:id + nx - 1), c(id:id + nx - 1), d(id:id + nx - 1), &
                       x(id:id + nx - 1), S(id:id + nx - 1))
        !print *, j
      end do

      do i = 1, nx
        do j = 1, ny
          sol(i, j) = sol(i, j) + rlx*x((j - 1)*nx + i)
          !sol(i,j) = x((j-1)*nx + i)
        end do
      end do

      do i = 1, nx
        do j = 1, ny

          id = (i - 1)*ny + j
          a(id) = w_d(i, j)
          b(id) = s_d(i, j)
          c(id) = o_d(i, j)
          d(id) = n_d(i, j)
          e(id) = e_d(i, j)

          S(id) = ind(i, j) - o_d(i, j)*sol(i, j)

          if (i > 1) S(id) = S(id) - w_d(i, j)*sol(i - 1, j)
          if (i < nx) S(id) = S(id) - e_d(i, j)*sol(i + 1, j)
          if (j > 1) S(id) = S(id) - s_d(i, j)*sol(i, j - 1)
          if (j < ny) S(id) = S(id) - n_d(i, j)*sol(i, j + 1)

        end do
      end do

      x = 0.0_DP
      ! do i = 1, nx
      !   do j = 1, ny
      !     x((i-1)*ny + j) = sol(i,j)
      !   end do
      ! end do

      do i = 1, nx

        id = (i - 1)*ny + 1

        if (i > 1) S(id:id + ny - 1) = S(id:id + ny - 1) - &
                                       a(id:id + ny - 1)*x(id - ny:id - 1)
        if (i < nx) S(id:id + ny - 1) = S(id:id + ny - 1) - &
                                        e(id:id + ny - 1)*x(id + ny:id + 2*ny - 1)
        call trisolver(b(id:id + ny - 1), c(id:id + ny - 1), d(id:id + ny - 1), &
                       x(id:id + ny - 1), S(id:id + ny - 1))

      end do

      do i = 1, nx
        do j = 1, ny
          sol(i, j) = sol(i, j) + rlx*x((i - 1)*ny + j)
          !sol(i,j) = x((i-1)*ny + j)
        end do
      end do
    end do
  end subroutine ADI

  subroutine gauss_seidel(w_d, e_d, s_d, n_d, o_d, ind, sol, error, swps, rlx)
    real(DP), intent(in) :: w_d(:, :), e_d(:, :), o_d(:, :), &
                            s_d(:, :), n_d(:, :), ind(:, :), rlx
    real(DP), intent(inout) :: sol(:, :), error
    integer, intent(in) :: swps
    real(DP), dimension(size(o_d, 1)*size(o_d, 2)) :: a, b, c, d, e, S, x

    integer :: nx, ny, i, j, k, id

    nx = size(o_d, 1)
    ny = size(o_d, 2)
    x = 0.0_DP
    do i = 1, nx
      do j = 1, ny

        id = (i - 1)*ny + j
        a(id) = s_d(i, j)
        b(id) = w_d(i, j)
        c(id) = o_d(i, j)
        d(id) = e_d(i, j)
        e(id) = n_d(i, j)

        S(id) = ind(i, j) - o_d(i, j)*sol(i, j)

        !x(id) = sol(i, j)
        if (j > 1) S(id) = S(id) - w_d(i, j)*sol(i, j - 1)
        if (j < ny) S(id) = S(id) - e_d(i, j)*sol(i, j + 1)
        if (i > 1) S(id) = S(id) - s_d(i, j)*sol(i - 1, j)
        if (i < nx) S(id) = S(id) - n_d(i, j)*sol(i + 1, j)
      end do
    end do
    do k = 1, swps

      !if (k == 1) print *, a
      !if (k == 1) print *, c
      do j = 1, ny
        do i = 1, nx
          id = (i - 1)*ny + j
          x(id) = S(id)
          if (j > 1) x(id) = x(id) - b(id)*x(id - 1)
          if (j < ny) x(id) = x(id) - d(id)*x(id + 1)
          if (i > 1) x(id) = x(id) - a(id)*x(id - ny)
          if (i < nx) x(id) = x(id) - e(id)*x(id + ny)
          x(id) = x(id)/c(id)
        end do
      end do
    end do
    error = norm2(x)
    do i = 1, nx
      do j = 1, ny
        id = (i - 1)*ny + j
        sol(i, j) = sol(i, j) + x(id)
        !sol(i,j) = x((j-1)*nx + i)
      end do
    end do
  end subroutine gauss_seidel
  subroutine trisolver(minusd, maind, plusd, hv, b)

    real(DP), intent(in) :: minusd(:), maind(:), &
                            plusd(:), b(:)
    real(DP), intent(inout) :: hv(:)
    real(DP) :: c(size(b)), d(size(plusd)), m, exact(size(hv))
    integer :: i, n
    n = size(hv)

    !print *, minusd(1), maind(1), plusd(1), b(1), &
    !minusd(n), maind(n), plusd(n), b(n)

    d(1) = plusd(1)/maind(1)
    c(1) = b(1)/maind(1)
    do i = 2, n
      !print *, minusd(i), maind(i), plusd(i), b(i)
      m = maind(i) - d(i - 1)*minusd(i)
      d(i) = plusd(i)/m
      c(i) = (b(i) - c(i - 1)*minusd(i))/m
    end do

    exact(n) = c(n)

    do i = n - 1, 1, -1
      exact(i) = c(i) - d(i)*exact(i + 1)
    end do

    hv = exact
  end subroutine trisolver

  subroutine write_tstep(field, param, t)
    type(point_t), intent(in) :: field(:, :)
    type(param_t), intent(in) :: param
    integer, intent(in) :: t

    integer :: i, j
    open (101, file="results.txt")!, position="append")

    do i = 1, param%nx
      do j = 1, param%ny
        if (j > param%my1 .and. j < param%my2) then
          if (i >= param%mx1 .and. i <= param%mx2) then
            write (101, "(6ES20.8, I3)") (t - 1)*param%ht, (i - 0.5_DP)*param%hx - param%xmax, &
              (j - 0.5_DP)*param%hy, field(i, j)%T, field(i, j)%F, field(i, j)%Z, 2
          else
            write (101, "(6ES20.8, I3)") (t - 1)*param%ht, (i - 0.5_DP)*param%hx - param%xmax, &
              (j - 0.5_DP)*param%hy, field(i, j)%T, field(i, j)%F, field(i, j)%Z, 3
          end if
        else
          write (101, "(6ES20.8, I3)") (t - 1)*param%ht, (i - 0.5_DP)*param%hx - param%xmax, &
            (j - 0.5_DP)*param%hy, field(i, j)%T, field(i, j)%F, field(i, j)%Z, 1
        end if

      end do
    end do
    write (101, *)
    write (101, *)
    close (101)
  end subroutine write_tstep
end module coefun_m
