module global_mod
    ! use mpi_f08 ! recommended over the mpif.h
    implicit none
    integer, parameter :: rp = selected_real_kind(15)
    real(kind = rp), parameter :: pi = 4.0_rp*atan(1.0_rp)
    ! for mpi:
    integer :: myRank, worldSize, ierror
    integer, dimension(:), allocatable :: workloads ! store workload for each core

end module global_mod

