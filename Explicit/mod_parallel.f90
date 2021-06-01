module mod_parallel
  
  ! Module that has parallel necesities

  use iso_fortran_env, only: int32, real32

  implicit none

  private

  public :: num_tiles, tiles_indices, tile_neighbors_1d, &
            tile_neighbors_2d

  interface tile_indices
    module procedure :: tile_indices_1d, tile_indices_2d
  end interface tile_indices
  contains

  pure function num_tiles(nprocs, nx, ny)
    ! Gives the number of images that are going to be
    ! either in vertical or horizontal

    integer(int32), intent(in   ) :: nprocs, nx, ny
    integer(int32), :: num_tiles(2)

    num_tiles(1) = max(int(sqrt(float(ny)/float(nx) * nprocs)), 1)

    num_tiles(2) = int(float(nprocs)/num_tiles(1))

    do while (num_tiles(1) * num_tiles(2) /= nprocs)
      num_tiles(1) = num_tiles(1) - 1
      num_tiles(2) = int(float(nprocs)/num_tiles(1))
    end do
  end function num_tiles

  pure function tile_indices_1d(dims, i, n) result(indices)
    ! Given global array size return start and end index

    integer(int32), intent(in   ) :: dims, i, n
    integer(int32) :: indices(2)
    integer(int32) :: offset, tile_size 

    tile_size = dims / n

    offset = mod(dims, n)
    
    if (i < offset) then 
      indices(1) = i * (tile_size + 1) + 1
      indices(2) = (i + 1) * (tile_size + 1)
    else 
      indices(1) = i * (tile_size + 1) + 1 + offset
      indices(2) = (i + 1) * (tile_size + 1) + offset
    end if

  end function tile_indices_1d

  pure function tiles_indices(dims) result(indices)
    ! Needs the number of points in both dimension [nx, ny]
    ! returns the start and finish indices for this images [i_s, i_e, j_s, j_e]

    integer(int32), intent(in   ) :: dims(2)
    integer(int32) :: indices(4)
    integer(int32) :: tiles(2), tiles_ij(2)

    tiles = num_tiles(num_images())
    tiles_ij = tilesn2ij(this_image())

    indices(1:2) = tiles_indices_1d(dims(1), tiles_ij(1), tiles(1))
    indices(3:4) = tiles_indices_1d(dims(2), tiles_ij(2), tiles(2))

  end function tiles_indices

  pure function tile_n2ij(n) result(ij)
    ! Give a index of the image and gets its position in the 2d layout
    ! j  _____________________________________________________________________
    !   |         6         |              7           |           8          |
    ! 2 |___________________|__________________________|______________________|
    !   |         3         |              4           |           5          |
    ! 1 |___________________|__________________________|______________________|
    !   |         0         |              1           |           2          |
    ! 0 |___________________|__________________________|______________________|
    !    i         0                       1                        2

    integer(int32), intent(in   ) :: n
    integer(int32) :: ij(2), i, j, tiles(2)

    if (n == 0) then
      ij = 0
    else
      tiles = num_tiles(num_images())
      i = mod(n, tiles(1))
      j = (n - 1) / tiles(1) 
      ij =  [i, j]
    end if

  end function tile_n2ij

  pure function tile_ij2n(ij) result(n)
    ! Give a position in the 2d layout and gets its index of the image 
    ! j  _____________________________________________________________________
    !   |         6         |              7           |           8          |
    ! 2 |___________________|__________________________|______________________|
    !   |         3         |              4           |           5          |
    ! 1 |___________________|__________________________|______________________|
    !   |         0         |              1           |           2          |
    ! 0 |___________________|__________________________|______________________|
    !    i         0                       1                        2

    integer(int32), intent(in   ) :: ij(2)
    integer(int32) :: n , i, j, tiles(2)

    tiles = num_tiles(num_images())

    n = (ij(2) - 1) * tiles(1) + ij(1)
  end function tile_ij2n

  pure function tile_neighbors_1d() result(neighbors)
    ! Returns indices of left and right neighbors

    integer(int32) :: neighbors(2)
    integer(int32) :: left, right

    if (num_images() > 1) then
      left = this_image() - 1
      right = this_image() + 1

      if (this_image() == 1) then
        left = num_images()
      else if (this_image() == num_images()) then
        right = 1
      end if
    else 
      left = 1
      right = 1      
    end if

    neighbors = [left, right]

  end function tile_neighbors_1d

  pure function tile_neighbors_2d(periodic) result(neighbors)
    ! Returns the neightbor image indices given
    logical, intent(in   ) :: preiodic
    integer(int32) :: neighbors(4)      
    integer(int32) :: tiles(2), tiles_ij(2), itile, jtile
    integer(int32) :: left, right, down, up
    integer(int32) :: ij_left(2), ij_right(2), ij_down(2), ij_up(2)

    tiles = num_tiles(num_images())
    !tiles_ij = tile_n2ij(this_image())
    [itile, jtile] = tile_n2ij(this_image())

    ij_left = [itile - 1, jtile]
    ij_right = [itile + 1, jtile]
    ij_down = [itile, jtile - 1]
    ij_up = [itile, jtile + 1]

    if (periodic) then

    else 

    end if 
  end function tile_neighbors_2d
end module mod_parallel