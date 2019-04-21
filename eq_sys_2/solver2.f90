module solve_methods
contains
  real, parameter :: eps = 10e-8

  !Функция, проверяющая диагональное преобладание в матрице
  logical function check_diagonal_dominance(A, n) result(res)
    real, dimension(:,:) :: A
    real :: tmp_sum
    logical ::  res
    integer :: n, i, j

    do i = 1,n
       tmp_sum = 0
       forall(j = 1:n, j .ne. i) tmp_sum = tmp_sum + abs(A(i,j))

       if(abs(a(i,i)) < tmp_sum) then
          res = .false.
          return
       end if
    end do
    res = .true.
    return
  end function check_diagonal_dominance

  !------------------Метод Якоби------------------------------!
  subroutine jakob_method(A, B, n, X)
    integer :: i, j, n
    real ::  norm_x ! норма разности X на к-м и (к-1)-м шаге
    real, dimension(:,:) :: A
    real, dimension(:) :: B, X
    real, allocatable, dimension(:) :: G
    real, allocatable, dimension(:,:) ::  Z

    if(check_diagonal_dominance(A,n) .ne. .true.) then
       write(*,*) "В матрице А нет диагонального преобладания"
    end if

    allocate(Z(n))
    allocate(G(n))


    ! Посчитаем матрицу Z = D^(-1)*(D − A)
    Z = 0
    forall(i = 1:n, j = 1:n, i .ne. j) Z(i,j) = -A(i,j)/A(j,j)

    !Посчитаем теперь ыектор G
    forall(i = 1:n) G(i) = B(i)/A(i,i)


    ! Приближаем пока norm_x = |X(k) − X(k−1)| > eps
    X = 0
    norm_x = 1
    do while(norm_x > eps)
       norm_x = sqrt(sum((matmul(X,Z)+G-X)**2))
       X = matmul(X,Z)+G
    end do

    deallocate(G)
    deallocate(Z)
  end subroutine jakob_method

  !------------------Метод Зейделя----------------------------!
  subroutine seidel_method(A,B,n,X)
    integer ::  i, j, n
    real ::  norm_x ! норма разности X на к-м и (к-1)-м шаге
    real, dimension(:,:) :: A
    real, dimension(:) :: B, X
    real, allocatable, dimension(:) :: Q, prevX
    real, allocatable, dimension(:,:) ::  P

    if(check_diagonal_dominance(A,n) .ne. .true.) then
       write(*,*) "В матрице А нет диагонального преобладания"
    end if
    allocate(Q(n))
    allocate(P(n,n))
    allocate(prevX(n))

    forall(i = 1:n, j = i:n) P(i,j) = -A(i,j)/A(i,i)
    forall(i = 1:n) Q(i) = B(i)/A(i,i)

    norm_x = 1
    X = 0
    do while(norm_x > eps)
       prevX = X
       do i = 1,n
          X(i) = Q(i) + sum(X(1:i-1)*P(i,1:i-1)) + sum(X(i+1:n)*P(i,i+1:n))
       end do

       norm_x = sqrt(sum((X-prevX)**2))
    end do
  end subroutine seidel_method

  !------------------Метод Релаксации----------------------------!
  subroutine relax_method(A,B,n,X)
    integer ::  i, j, n
    real ::  maxcoord
    real, dimension(:,:) :: A
    real, dimension(:) :: B, X
    real, allocatable, dimension(:) :: Q
    real, allocatable, dimension(:,:) ::  P

    if(check_diagonal_dominance(A,n) .ne. .true.) then
       write(*,*) "В матрице А нет диагонального преобладания"
    end if

    allocate(Q(n))
    allocate(P(n,n))

    forall(i = 1:n, j = i:n) P(i,j) = -A(i,j)/A(i,i)
    forall(i = 1:n) Q(i) = B(i)/A(i,i)

    X = 0 ! Начальое приближение !

    do while(maxval(abs(Q)) > eps)
       maxcoord = maxloc(abs(Q))
       X(maxcoord) = x(maxcoord) + Q(maxcoord)

       ! Вычисляем новые невязки !
       forall(i = 1:n) Q(i) = Q(i) + P(i,maxcoord)*Q(maxcoord)
    end do
  end subroutine relax_method
