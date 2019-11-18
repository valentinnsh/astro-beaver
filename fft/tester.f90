program tester_fft
  use fft
  implicit none
  
  integer :: i, N, N_2
  character(3) :: key
  complex, allocatable, dimension(:) :: X, Y, W
  real, allocatable, dimension(:, :) :: data_X, data_Y  
  
  open(10, file = 'data.dat')
  read(10,'(2x,I6)') N

  N_2 = 2**ceiling(log(real(N))/log(2.0))

  allocate(data_X(0:N_2-1, 2))

  do i = 0, N - 1
     read(10,*) data_X(i, :)
  end do
     
  do i = N, N_2 - 1
     data_X(i, :) = 0
  end do
  
  N = N_2

  close(10)

  allocate(X(0:N-1), Y(0:N-1), data_Y(0:N-1, 2), W(0:N-1))

  ! Делаем из двух столбцов массива data_X один комплекснозначный массив X

  X = cmplx(data_X(:, 1), data_X(:, 2))

  deallocate(data_X)

  call getarg(1, key)
  select case (key)
  case('dir')
     call  W_negative(W, N)
     call bit_reversal_sorting(X,Y,W)
  case('inv')
     call W_positive(W, N)
     call bit_reversal_sorting(X,Y,W)
  case default
     stop'Error, wrong key'
  end select

  Y = Y/sqrt(real(N)) ! Нормируем Y
  deallocate(X)

  data_Y(:, 1) = real(Y)
  data_Y(:, 2) = aimag(Y)
  
  open(20, file = 'result.dat')
  write(20, *) '# ', N
  do i = 0, N - 1
     write(20, *) data_Y(i, :)
  end do
  close(20)
  
  open(30, file = 'abs.dat')
  do i = 0, N - 1
     write(30, *) abs(Y(i))
  end do
  close(30)

  deallocate(Y, data_Y)  
  
end program tester_fft
