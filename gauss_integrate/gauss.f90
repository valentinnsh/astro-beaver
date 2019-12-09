module gauss
  use legendre
  use eq_solver
  use file_ops
  implicit none
contains

  function gauss_quad_coef(n, t) result(A)
    real(8) :: t(1:n)
    integer(4) :: n
    real(8), allocatable :: A(:), M(:,:), B(:)
    integer(4) :: i, k
    character :: key = 'i'

    allocate(A(1:n), M(1:n, 1:n), B(1:n))

    call quicksort(t,1,n)

    do k = 1, n
       do i = 1, n
          M(k, i) = t(i)**(k - 1.0_8)
       end do
    end do

    B = 0.0_8

    do k = 1, n, 2
       B(k) = 2.0_8/k
    end do

    call find_solution(M, B, key, n, A)
  end function  gauss_quad_coef

  function integration(lim_a, lim_b, n, f) result(res)
    real(8) :: lim_a, lim_b, res
    character(len=14) :: filename
    integer(4) :: i, n
    real(8), dimension(1:n) :: A, t

    interface
       function f(x)
         real(8) :: x, f
       end function f
    end interface

    filename = get_file_name(n)

    open(10, file=filename)
    do i=1,n
       read(10,*) A(i), t(i)
    enddo
    close(10)

    res = 0

    do i = 1, n
       res = res + A(i)*(lim_b-lim_a)/2.0_8 * f( (t(i) * (lim_b - lim_a) / 2.0_8 + (lim_a + lim_b) / 2.0_8 ))
    end do
  end function integration

end module gauss
