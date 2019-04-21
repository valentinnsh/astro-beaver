program testing_prog
  use solve_methods
  use matrixop
  implicit none

  integer :: i, j, n
  real, allocatable, dimension(:,:) :: A
  real, allocatable, dimension(:) :: B, X, R
  real :: R_module
  ! Вводится с комндной строки
  ! jakobe, seidel или relax
  character*6 :: method
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

  allocate(X(n))

  select case(method)
  case("jakobe")
     call jakob_method(A,B,n,X)
  case("seidel")
     call seidel_method(A,B,n,X)
  case("relax ")
     call relax_method(A,B,n,X)
  case default
     write(*,*) "Error, wrong method name.\n"
  end select

  open(200, file = 'result.dat')
  write(200,*) "# ", n
  do i = 1,n
     write(200,*) X(i)
  end do

  ! Ищем вектор невязки
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
end program testing_prog
