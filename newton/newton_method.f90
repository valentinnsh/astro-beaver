
module newton_method_func
  use eq_solver

  implicit none

contains

  ! пример функции которую будем передавать в программу !
  function myfun(xin) result(xout)
    real, dimension(:) :: xin
    real, dimension(1:size(xin)) :: xout
    integer :: i

    do i = 1, size(xin)
       xout(i) = xin(i)**i - 14**i
    end do
    return
  end function myfun



  ! вычисление якобиана !
  subroutine prive_diff(x_orig, df, func_name)
    external func_name
    interface
       function func_name(xin) result(xout)
         real, dimension(:) :: xin
         real, dimension(1:size(xin)) :: xout
         integer :: i
       end function func_name
    end interface

    real, dimension(:) :: x_orig
    real, dimension(1:size(x_orig)) :: x
    real, dimension(1:size(x_orig),1:size(x_orig)), intent(out) :: df
    integer :: i, n
    real :: eps = 1e-3

    n = size(x_orig)
    do i = 1, n
       x = x_orig
       x(i) = x_orig(i) + eps
       df(i,1:n) = (func_name(x) - func_name(x_orig))/eps
    end do

  end subroutine prive_diff



  function newton_method(fname, xin, max_iter) result(res)
    external fname
    interface
       function fname(xin) result(xout)
         real, dimension(:) :: xin
         real, dimension(1:size(xin)) :: xout
         integer :: i
       end function fname
    end interface

    integer :: i, max_iter, n
    real, dimension(:) :: xin
    real, dimension(1:size(xin)) ::  x, xprev, B, res
    real, dimension(1:size(xin),1:size(xin)) :: A, df
    real :: eps = 1e-4
    xprev = xin
    x = xprev + 2.1
    n = size(xin)

    i = 1
    do while((sum(abs(x-xprev)) > eps) .and. i < max_iter)

       xprev = x
       call prive_diff(xprev,df,fname)
       A(1:n, 1:n) = -df(1:n,1:n)
       B = fname(xprev)
       call find_solution(A, B, 'l', n, x)
       x = x + xprev
       i = i + 1
    end do

    res = x
  end function newton_method


end module newton_method_func
