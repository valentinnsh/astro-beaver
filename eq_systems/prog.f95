program solve_equalation
  use eq_solver
  implicit none

  integer :: id, i, j, n
  real, allocatable, dimension(:,:) :: A
  real, allocatable, dimension(:) :: B, X, R
  real :: R_module
  !Вводится с комндной строки
  ! g для схемы гаусса, j для схемы жордана,
  ! l для схемы гаусса с перестановкой строк
  character :: method
  !получаем название метода, с помощью которого будем находить решениец
  call getarg(1,method)

  open(100, file = 'data.dat')
  read(100,'(2x,I6)') n
  allocate(A(n,n))
  !Считываем матрицу коэффициентов
  call read_matrix(100, n, A)

  allocate(B(n))
  !Считаем вектор свободных значений
  do i=1,n
     read(100,*) B(i)
  end do

  close(100)


  !X - вектор решений
  allocate(X(n))
  X = 0

  call find_solution(A,B,method,n,X)


  !Вывод результата в соответствующий файл
  open(200, file = 'result.dat')
  write(200,*) "# ", n
  do i = 1,n
     write(200,*) X(i)
  end do

  !Ищем вектор невязки
  allocate(R(n))
  R = 0
  do i = 1,n
     do j = 1,n
        R(i) = R(i) + A(i,j)*X(j)
     end do
     R(i) = R(i) - B(i)
  end do
  R_module = 0
  do i = 1,n
     R_module = R_module + R(i)**2
  end do

  R_module = sqrt(R_module)
  write(*,*) R_module
end program solve_equalation
