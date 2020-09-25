module paralel_m

  use type_m, only: DP

  implicit none 
  
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
    contains

    procedure :: read_parll
  end type parll_t
  contains
  
  subroutine read_parll(this)

    class(parll_t) :: this

    !open(101, file = "dataparll", status = "old")
    !  read(101,*) this%nrows, this%ncols
    !close(101)
    end subroutine
end module paralel_m