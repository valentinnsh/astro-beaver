program tester_fft
  use fft

  integer :: i, N;
  complex, allocatable, dimension(:,:) :: x, y
  open(10, file = 'data.dat')
  open(20, file = 'result.dat')

  read(10,'(2x,I6)') N
  allocate(x(0:N-1, 2)); allocate(y(0:N-1, 2))

  do i = 1, n+1
     read(10,*) XYP(i,1:3)
  end do

end program tester_fft
