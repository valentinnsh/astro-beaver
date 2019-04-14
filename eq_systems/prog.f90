program solve_equalation
  use eq_solver
  implicit none

  integer :: id, i, n
  real, dimension(:,:) :: A
  real, dimension(:) :: B, X

  !Вводится с комндной строки
  ! gaus для схемы гаусса, jordan для схемы жордана,
  ! upgaus для схемы гаусса с перестановкой строк
  character(len = 6) :: method
  !получаем название метода, с помощью которого будем находить решениец
  call getarg(1,method)

  open(100, file = 'data.dat')
  read(100,'(2x,I6)') n
  allocate(A(n,n))
  !Считываем матрицу коэффициентов
  call read_matrix(100, n, A)

  allocate(B(n))
  !Считаем вектор свободных значений
  forall(i=1:n) read(100,*) B(i)
  close(100)

  !X - вектор решений
  allocate(X(n))
  X = 0

  call find_solution(A,B,method,n,X)

  !Вывод результата в соответствующий файл
  open(200, file = 'result.dat')
  write(200,*) "# ", n
  forall(i = 1:n) write(100,*) X(i)

end program solve_equalation
