program triple_diag_matrix
  use tri_makedata
  use tri_mult
  implicit none
  real, allocatable, dimension(:,:) :: matrix_A, matrix_B, matrix_C
  integer :: n, i,j

  !call make_data(10)
  
  open(100, file='data1.dat')
  read(100,'(2x,I6)') n
  allocate(matrix_A(n,3))
  read(100,*)matrix_A(1,2:3)
  matrix_A(1,1)=0
  do i=2,n-1
     read(100,*)matrix_A(i,1:3)
  enddo
  read(100,*)matrix_A(n,1:2)
  matrix_A(n,3)=0
  close(100)
  
  open(100, file='data2.dat')
  read(100,'(2x,I6)') n
  allocate(matrix_B(n,3))
  read(100,*)matrix_B(1,2:3)
  matrix_B(1,1)=0
  do i=2,n-1
     read(100,*)matrix_B(i,1:3)
  enddo
  read(100,*)matrix_B(n,1:2)
  matrix_B(n,3)=0
  close(100)
  
  allocate(matrix_C(n,5))
  matrix_C=0
  
  call tri_matrix_mult(n, matrix_A, matrix_B, matrix_C)
  
  open(100, file='result.dat')
  write(100,"('# ',I6)") n
  write(100,*)matrix_C(1,3:5)
  write(100,*)matrix_C(2,2:5)
  do i=3,n-2
     write(100,*)matrix_C(i,1:5)
  enddo
  write(100,*)matrix_C(n-1,1:4)
  write(100,*)matrix_C(n,1:3)   
  close(100)
    
  deallocate(matrix_A)
  deallocate(matrix_B)
  deallocate(matrix_C)

  
end program triple_diag_matrix

