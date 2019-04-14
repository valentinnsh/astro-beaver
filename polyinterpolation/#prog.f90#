program poly_interpolation
  use interpolation
  implicit none
  
!Считываем тип сетки, для которой нужно вычислить значение 
character(len=9) :: net
integer n, i
real a,b
! введем две матрицы. В N+1x2 матрице func_data будут находиться значения интерполируемой
! функции и вычесленные в соответствии с выбранной сеткой значения х
! А в матрице 100N+1х2 inter_poly будем класть результат - значения интерполяционного полинома
! на равномерной сетке
real, allocatable, dimension(:,:) :: func_data, poly_data

call getarg(1,net)

if (net == 'chebyshev') then
   open(100, file = 'chebyshev.dat')
   read(100,'(2x,I6)') n
   read(100,*) a, b
   allocate(func_data(0:n,1:2))
   func_data = 0
   allocate(poly_data(0:100*n,2))
   poly_data = 0
   read(100,*) func_data(0:n,2)
   close(100)
   
   open(200, file = 'res_chebyshev.dat')
   call chebyshev_net(a,b,n,func_data)
   call interpolate(a,b,n,func_data,poly_data)

else
   open(100, file = 'uniform.dat')
   read(100,'(2x,I6)') n
   read(100,*) a, b
   allocate(func_data(0:n,1:2))
   func_data = 0
   allocate(poly_data(0:100*n,2))
   poly_data = 0
   read(100,*) func_data(0:n,2)
   close(100)
   
   open(200, file = 'res_uniform.dat')
   call uniform_net(a,b,n,func_data)
   call interpolate(a,b,n,func_data,poly_data)

end if

do i = 0, 100*n+1
   write(200,*) poly_data(i,1:2)
end do

close(200)

deallocate(func_data)
deallocate(poly_data)

end program
