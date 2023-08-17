      program mpi_ver
      implicit none
      include 'mpif.h'
      integer(kind=kind(MPI_VERSION)), parameter :: zero = ichar('0')
      character, dimension(17), parameter :: mpiver_str =&
      (/ 'I', 'N', 'F', 'O', ':', 'M', 'P', 'I', '-', 'V', 'E', 'R', '[', &
        char(zero + MPI_VERSION), &
        '.', &
        char(zero + MPI_SUBVERSION), ']' /)
      print *, mpiver_str
      end program mpi_ver
