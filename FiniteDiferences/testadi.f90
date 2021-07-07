program test_ADI
  use type_m, only: DP
  use coefun_m, only: ADI, gauss_seidel
  implicit none

  include 'mpif.h'

  ! THERE IS 12 POINTS
  ! the mesh is 4x3
  !
  !  8  9 10 11
  !  4  5  6  7
  !  0  1  2  3
  !
  !   3   1   0   0   2   0   0   0   0   0   0   0  !!!!!!   0.5
  !   1   3   1   0   0 0.5   0   0   0   0   0   0  !!!!!!   0.5
  !   0   1   2   2   0   0 0.6   0   0   0   0   0  !!!!!!   0.5
  !   0   0   2   5   0   0   0 0.7   0   0   0   0  !!!!!!   0.5
  !   2   0   0   0   4   1   0   0 0.8   0   0   0  !!!!!!   0.5
  !   0 0.5   0   0   1   6   1   0   0 0.9   0   0  !!!!!!   0.5
  !   0   0 0.6   0   0   1   3   1   0   0 0.9   0  !!!!!!   0.5
  !   0   0   0 0.7   0   0   1   2   0   0   0   1  !!!!!!   0.5
  !   0   0   0   0 0.8   0   0   0   7   2   0   0  !!!!!!   0.5
  !   0   0   0   0   0 0.9   0   0   2   8   2   0  !!!!!!   0.5
  !   0   0   0   0   0   0 0.9   0   0   2   4   1  !!!!!!   0.5
  !   0   0   0   0   0   0   0   1   0   0   1   3  !!!!!!   0.5

  real(DP), dimension(:, :), allocatable :: w, e, o, s, n, ind, sol
  real(DP) :: rlx = 1.0, error
  real(DP) :: A(12, 12), b(12), realsol(12)
  real(dp), dimension(12) :: realtrisol, trisol, md, maind, pd

  integer :: ierr, nprocs, myproc, i, maxit = 10000, swps = 1
  integer :: status(MPI_STATUS_SIZE)
  integer, parameter :: TOP_TAG = 13, BOT_TAG = 14

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myproc, ierr)

  A(1, :) = (/3.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  A(2, :) = (/1.0, 3.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  A(3, :) = (/0.0, 1.0, 2.0, 2.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0/)
  A(4, :) = (/0.0, 0.0, 2.0, 5.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0/)
  A(5, :) = (/2.0, 0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
  A(6, :) = (/0.0, 0.5, 0.0, 0.0, 1.0, 6.0, 1.0, 0.0, 0.0, 0.9, 0.0, 0.0/)
  A(7, :) = (/0.0, 0.0, 0.6, 0.0, 0.0, 1.0, 3.0, 1.0, 0.0, 0.0, 0.9, 0.0/)
  A(8, :) = (/0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0/)
  A(9, :) = (/0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 7.0, 2.0, 0.0, 0.0/)
  A(10, :) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 2.0, 8.0, 2.0, 0.0/)
  A(11, :) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 2.0, 4.0, 1.0/)
  A(12, :) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 3.0/)
  b = 0.5_dp
  realsol = 0.0_DP

  if (nprocs == 3) then
    allocate (w(3, 4), e(3, 4), o(3, 4), s(3, 4), n(3, 4), ind(3, 4), sol(3, 4))
    w(:, :) = 0.0
    e(:, :) = 0.0
    s(:, :) = 0.0
    n(:, :) = 0.0
    o(:, :) = 1.0
    sol(:, :) = 0.0
    ind(:, :) = 0.0
    if (myproc == 0) then
      w(2, :) = (/0.0, 1.0, 1.0, 2.0/)
      e(2, :) = (/1.0, 1.0, 2.0, 0.0/)
      o(2, :) = (/3.0, 3.0, 2.0, 5.0/)
      s(2, :) = (/0.0, 0.0, 0.0, 0.0/)
      n(2, :) = (/2.0, 0.5, 0.6, 0.7/)
    else if (myproc == 1) then
      w(2, :) = (/0.0, 1.0, 1.0, 1.0/)
      e(2, :) = (/1.0, 1.0, 1.0, 0.0/)
      o(2, :) = (/4.0, 6.0, 3.0, 2.0/)
      s(2, :) = (/2.0, 0.5, 0.6, 0.7/)
      n(2, :) = (/0.8, 0.9, 0.9, 1.0/)
    else if (myproc == 2) then
      w(2, :) = (/0.0, 2.0, 2.0, 1.0/)
      e(2, :) = (/2.0, 2.0, 1.0, 0.0/)
      o(2, :) = (/7.0, 8.0, 4.0, 3.0/)
      s(2, :) = (/0.8, 0.9, 0.9, 1.0/)
      n(2, :) = (/0.0, 0.0, 0.0, 0.0/)
    end if
    ind(2, :) = (/0.5, 0.5, 0.5, 0.5/)
  else if (nprocs == 4) then
    allocate (w(3, 1), e(3, 1), o(3, 1), s(3, 1), n(3, 1), ind(3, 1), sol(3, 1))
  else if (nprocs == 1) then
    allocate (w(3, 4), e(3, 4), o(3, 4), s(3, 4), n(3, 4), ind(3, 4), sol(3, 4))
    w(1, :) = (/0.0, 1.0, 1.0, 2.0/)
    w(2, :) = (/0.0, 1.0, 1.0, 1.0/)
    w(3, :) = (/0.0, 2.0, 2.0, 1.0/)
    e(1, :) = (/1.0, 1.0, 2.0, 0.0/)
    e(2, :) = (/1.0, 1.0, 1.0, 0.0/)
    e(3, :) = (/2.0, 2.0, 1.0, 0.0/)
    o(1, :) = (/3.0, 3.0, 2.0, 5.0/)
    o(2, :) = (/4.0, 6.0, 3.0, 2.0/)
    o(3, :) = (/7.0, 8.0, 4.0, 3.0/)
    s(1, :) = (/0.0, 0.0, 0.0, 0.0/)
    s(2, :) = (/2.0, 0.5, 0.6, 0.7/)
    s(3, :) = (/0.8, 0.9, 0.9, 1.0/)
    n(1, :) = (/2.0, 0.5, 0.6, 0.7/)
    n(2, :) = (/0.8, 0.9, 0.9, 1.0/)
    n(3, :) = (/0.0, 0.0, 0.0, 0.0/)
    ind(:, :) = 0.5
    sol(:, :) = 0.0
  else
    print *, "The number of processor to test has to be either 1, 3 or 4."
    stop
  end if

  call iterative_gauss_seidel(a, b, realsol, 10d-8, maxit, i)

  trisol = 0.0_DP
  realtrisol = 0.0_DP
  md = 0.0_DP
  pd = 0.0_DP

  do i = 1, 12
    maind(i) = A(i, i)
    if (i > 1) md(i) = A(i, i - 1)
    if (i > 2) A(i, 1:i - 2) = 0.0_DP
    if (i < 10) A(i, i + 2:12) = 0.0_DP
    if (i < 12) pd(i) = A(i, i + 1)
  end do
  print *, "----------------------------------------------"
  call iterative_gauss_seidel(a, b, realtrisol, 10d-8, maxit, i)
  call trisolver(md, maind, pd, trisol, b)
  print *, abs(trisol - realtrisol)
  maxit = 10

  call adi(w, e, s, n, o, ind, sol, error, 4, rlx)
  do i = 1, 10000! maxit
    if (nprocs == 3) then
      call gauss_seidel(w, e, s, n, o, ind, sol, error, swps, rlx)
      if (myproc == 0) then
        call MPI_SENDRECV(sol(2, :), size(sol(2, :)), MPI_DOUBLE_PRECISION, &
                          1, TOP_TAG, sol(3, :), size(sol(3, :)), MPI_DOUBLE_PRECISION, &
                          1, TOP_TAG, MPI_COMM_WORLD, status, ierr)
        ind(3, :) = sol(3, :)
      elseif (myproc == 1) then
        call MPI_SENDRECV(sol(2, :), size(sol(2, :)), MPI_DOUBLE_PRECISION, &
                          0, TOP_TAG, sol(1, :), size(sol(1, :)), MPI_DOUBLE_PRECISION, &
                          0, TOP_TAG, MPI_COMM_WORLD, status, ierr)
        call MPI_SENDRECV(sol(2, :), size(sol(2, :)), MPI_DOUBLE_PRECISION, &
                          2, TOP_TAG, sol(3, :), size(sol(3, :)), MPI_DOUBLE_PRECISION, &
                          2, TOP_TAG, MPI_COMM_WORLD, status, ierr)
        ind(1, :) = sol(1, :)
        ind(3, :) = sol(3, :)
      elseif (myproc == 2) then
        call MPI_SENDRECV(sol(2, :), size(sol(2, :)), MPI_DOUBLE_PRECISION, &
                          1, TOP_TAG, sol(1, :), size(sol(1, :)), MPI_DOUBLE_PRECISION, &
                          1, TOP_TAG, MPI_COMM_WORLD, status, ierr)
        ind(1, :) = sol(1, :)
      end if
    else if (nprocs == 1) then
      call gauss_seidel(w, e, s, n, o, ind, sol, error, swps, rlx)
    end if

    !if (mod(i, maxit/10) == 0) print *, myproc, ":", error
  end do
  !print *, realsol
  if (nprocs == 3) then
    if (myproc == 0) then
      print *, myproc, ":", abs(sol(2, :) - realsol(1:4))!(/0.05274725274725275, -0.013651877133105802, 0.4597315436241611, -0.18900343642611683/))
    else if (myproc == 1) then
      print *, myproc, ":", abs(sol(2, :) - realsol(5:8))!(/0.17771084337349397, 0.05701754385964912, -0.04633204633204633, 0.243006993006993/))
    else if (myproc == 2) then
      print *, myproc, ":", abs(sol(2, :) - realsol(9:12))!(/0.005813953488372093, 0.037091988130563795, 0.07017543859649122, 0.18681318681318682/))
    end if
  else if (nprocs == 1) then
    print *, "----------------------------------------------------------"
    print *, myproc, ":", abs(sol(1, :) - realsol(1:4))!(/0.19506172839506172, -0.056856187290969896, 0.4217252396166134, -0.11312217194570136/))
    print *, myproc, ":", abs(sol(2, :) - realsol(5:8))!(/-0.01417004048582996, 0.10756972111553785, -0.10059171597633136, 0.3173431734317343/))
    print *, myproc, ":", abs(sol(3, :) - realsol(9:12))!(/0.07374631268436578, -0.0024330900243309003, 0.1375770020533881, 0.045081967213114756/))
  end if
  call MPI_FINALIZE(ierr)
contains
  SUBROUTINE iterative_gauss_seidel(a, b, x, tol, iter, clave)
    USE, INTRINSIC :: iso_fortran_env, ONLY: WP => REAL64
    IMPLICIT NONE

    ! -----------------------------------------------------------------
    ! Argumentos de la subrutina
    ! -----------------------------------------------------------------
    REAL(WP), INTENT(IN) :: a(:, :) ! Matriz del sistema
    REAL(WP), INTENT(IN) :: b(:) ! Término independiente
    REAL(WP), INTENT(INOUT) :: x(:) ! Aproximación inicial/solución
    REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error
    INTEGER, INTENT(INOUT) :: iter ! Max iteraciones/iter realizadas
    INTEGER, INTENT(OUT) :: clave ! Clave de éxito/error
    ! clave = 0, OK.
    ! clave /= 0, max iter. alcanzado
    ! -----------------------------------------------------------------
    ! Variables locales
    ! -----------------------------------------------------------------
    INTEGER :: ndim ! Dimensión del problema
    INTEGER :: k ! Contador de iteraciones
    INTEGER :: j ! Indice
    REAL(WP) :: xi ! Componente del vector iterado
    REAL(WP) :: xnorma ! Norma del vector iterado
    REAL(WP) :: difnorma ! Norma de la diferencia entre dos
    ! iteraciones del vector
    ! -----------------------------------------------------------------
    ! Procedimiento
    ! -----------------------------------------------------------------
    clave = 1
    ndim = SIZE(a, 1)
    DO k = 1, iter
      xnorma = 0.0_WP
      difnorma = 0.0_WP
      DO j = 1, ndim
        xi = (b(j) - DOT_PRODUCT(a(j, 1:j - 1), x(1:j - 1)) &
        & - DOT_PRODUCT(a(j, j + 1:ndim), x(j + 1:ndim)))/a(j, j)
        xnorma = MAX(xnorma, ABS(xi))
        difnorma = MAX(difnorma, ABS(xi - x(j)))
        x(j) = xi
      END DO
      IF (difnorma <= xnorma*tol) THEN
        iter = k
        clave = 0
        EXIT
      ENDIF
    ENDDO
  END SUBROUTINE iterative_gauss_seidel

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
end program test_ADI
