
program newton_method_tester
  use newton_method_func
  implicit none

  integer :: i, j, n, max_iter
  real, allocatable, dimension(:) :: xin, solvation
  real :: abs_x
  !   Возьмем к примеру систему из 6 уравнений   !
  ! А макс. колличество итераций пусть будет 500 !

  n = 3; max_iter = 500
  allocate(xin(n)); allocate(solvation(n))
  abs_x = 0
  xin = 1.0

  solvation = newton_method(myfun, xin, max_iter)
  do i = 1,n
     abs_x = abs_x + solvation(i)**2
  end do
  abs_x = abs_x**0.5

  write(*,*) abs_x

  open(20, file = 'result.dat')
  do i = 1,n
     write(20, *) solvation(i)
  end do
  close(20)
end program NEWTON_METHOD_TESTER
