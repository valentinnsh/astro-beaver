program task1
  use matrix
  use matrix_dat
  implicit none
  real, allocatable, dimension(:,:) :: matrix_A, matrix_B, matrix_C
  integer :: n
  real :: start, finish

  call make_data(10)
  
  open(1, file='data1.dat')
  read(1,'(2x,I6)') n
  allocate(matrix_A(n,n))

  call read_matrix(1, n, matrix_A)
  close(1)
  
  open(2, file='data2.dat')
  allocate(matrix_B(n,n))
  read(2,'(2x,I6)') n

  call read_matrix(2,n,matrix_B)
  close(2)

  allocate(matrix_C(n,n))

  call native_mult(n, matrix_A, matrix_B, matrix_C)

  open(10, file='result.dat')
  write( 10,"('# ',I6)") n
  call print_matrix(10, n, matrix_C)
  close(10)

  deallocate(matrix_A)
  deallocate(matrix_B)
  deallocate(matrix_C)

end program task1

