module mod_common_mpi
  use mpi
  use mod_param, only: dims
  implicit none
  integer :: myid
  integer, dimension(0:1,3) :: nb
  logical, dimension(0:1,3) :: is_bound
  integer, dimension(3) :: coord
  integer :: comm_cart,ipencil,ierr
  integer, dimension(3) :: halo
  integer :: status(MPI_STATUS_SIZE)
  integer, dimension(3,3) :: dims_xyz
  integer, dimension(3) :: ijk_min
end module mod_common_mpi
