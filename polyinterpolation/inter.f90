module interpolation
contains
  subroutine uniform_net(a,b,n,func_data)
    integer n, i
    real a,b
    real, allocatable, dimension(:,:) :: func_data
    
    do i = 0, n
       func_data(i,1) = a + (b-a)/n
    end do
  end subroutine uniform_net

  subroutine chebyshev_net(a,b,n,func_data)
    integer n, i
    real a,b, pi
    real, allocatable, dimension(:,:) :: func_data

    pi = 4*atan(1.0)
    do i = 0, n
       func_data(i,1) = (a+b)/2 + (a-b)/2*cos((2*i+1)*pi/(2*n+2))
    end do
  end subroutine chebyshev_net
  
  
  subroutine interpolate(a,b,n,func_data,poly_data)
    integer n, n1, i
    real a,b, dx, curr_x
    real, allocatable, dimension(:,:) :: func_data, poly_data

    n1 = n*100
    dx = (b-a)/n1
    curr_x = a
    
    do i = 0, n1
       poly_data(i,1) = curr_x
       poly_data(i,2) = calc_poly(curr_x, n, func_data)
       curr_x = curr_x + dx
    enddo
  end subroutine interpolate

  !функция считает значение полинома Лагранжа для функции заданной в func_data в конкретной точке
  function calc_poly(curr_x, n, func_data)  result(res)
    integer :: n,i, k
    real curr_x, res, base_k
    real, allocatable, dimension(:,:) :: func_data
    res = 0
    
    do k = 0, n
       base_k = 1
       do i = 0, n
          if(i .ne. k) then
             base_k = base_k*(curr_x-func_data(i,1))/(func_data(k,1) - func_data(i,1))
          end if
       end do
       res = res+func_data(k,2)*base_k
    end do
  end function calc_poly
  
end module interpolation
