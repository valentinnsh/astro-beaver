module horner_and_bernoulli
  implicit none
contains

  subroutine eq_solver_bernoulli(A0, X)
    real(8), dimension(0:) :: A0
    real(8), dimension(0:(size(A0) - 1)/2) :: A, B
    real(8), dimension(1:size(A0)-1) :: X
    real(8), dimension(1:size(A0)-1) :: Y
    real(8) :: root, sqroot
    integer(4) :: n, i

    n = size(A0) - 1

    if (mod(n,2) .ne. 0) then
        X(n) = 0.0_8
    endif

    A(0:n/2) = A0(0:n:2) ! Т.к. в полиноме Лежандра есть либо чётные степени, либо нечётные


    do i = 0, n/2 - 3
        call RANDOM_NUMBER(Y)
        call find_root(A(0:n/2-i), Y(1:n/2 - i), root)

        X(2*i + 1) = sqrt(root)
        X(2*i + 2) = -sqrt(root)

        call horner_division(A(0:n/2-i), B(0:n/2-i), root)
        A(0:n/2-i-1) = B(0:n/2-i-1)
    enddo

    sqroot = sqrt(A(1)**2.0_8 - 4.0_8*A(0)*A(2))
    X(2*i + 1) = sqrt((-A(1) + sqroot)/(2.0_8*A(0)))
    X(2*i+2) = - X(2*i+1)
    X(2*i+3) = sqrt((-A(1) - sqroot)/(2.0_8*A(0)))
    X(2*i+4) = - X(2*i+3)
  end subroutine eq_solver_bernoulli

  ! Деление многочлена на многочлен методом Горнера
  subroutine horner_division(num, den, root)
    real(8), dimension(0:) :: num
    real(8) :: root
    real(8), dimension(0:size(num) - 1) :: den
    integer(4) :: i, n

    n = size(num) - 1
    den(0) = num(0)

    do i = 1, n
       den(i) = den(i - 1)*root + num(i)
    end do

  end subroutine horner_division

  ! Ищем максимальный корень для метода Бернулли
  subroutine find_root(A, Y, root)
    real(8), dimension(0:) :: A
    real(8) :: root
    real(8), dimension(1:size(A)) :: Y
    integer(4) :: i, n
    real(8) :: eps = 0.000000001_8

    n = size(A) - 1

    do while (abs(Y(n)/Y(n - 1) - Y(n - 1)/Y(n - 2)) > eps)
       Y(n+1) = 1.0_8/(-A(0))*sum(A(1:n)*Y(n:1:-1)) ! Итеративная формула для вычисления Y(n+1)
       forall (i = 1:n-1) Y(i) = Y(i+1)
       Y(n) = Y(n+1)
    end do

    root = Y(n)/Y(n - 1)

  end subroutine find_root


end module horner_and_bernoulli
