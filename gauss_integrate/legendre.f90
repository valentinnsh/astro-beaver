module legendre
  use horner_and_bernoulli
    implicit none

    contains

    ! Вычисляем коэффициенты Лежандра
    recursive function coeff_legendre(n) result(P)
    integer(4) :: n
    real(8), dimension(0:n) :: P
    real(8), dimension(0:n-1) :: P1
    real(8), dimension(0:n-2) :: P2
    integer(4) :: i


    if (n==0) then
        P(0) = 1.0_8
    else if (n==1) then
        P(0:1) = (/ 1.0_8, 0.0_8 /)
    else
        P1 = coeff_legendre(n-1)
        P2 = coeff_legendre(n-2)

        P(0) = (2.0_8*n-1.0_8)/n * P1(0)
        P(1) = (2.0_8*n-1.0_8)/n * P1(1)

        forall (i=2:n-1) P(i) = (2.0_8*n-1.0_8)/n*P1(i)-(n-1.0_8)/n*P2(i-2)
        P(n) = - (n - 1.0_8)/n*P2(n-2)
    end if
  end function coeff_legendre


    ! Ищем  корни полиномов Лежандра
    function legendre_roots(n) result(X)
    integer(4) :: n
    real(8), dimension(1:n) :: X
    real(8), dimension(0:n) :: A

    select case(n)
        case(1); X=0.0_8 ! P1 = x
        case(2); X=(/ 0.57735026919_8, -0.57735026919_8 /) ! P2 = (3*x^2 -1)/2
        case(3); X=(/ 0.77459666924_8, 0.0_8, -0.77459666924_8 /) ! P3 = (5*x^3 - 3*x)/2
        case default
            A =  coeff_legendre(n)
            call eq_solver_bernoulli(A, X)
    end select

  end function legendre_roots
end module legendre
