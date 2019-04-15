module eq_solver
contains
  subroutine print_matrix(id, n, matrix)
    integer :: id, i, n
    real, dimension(:,:) :: matrix
    do i = 1, n
       write(id, *) matrix(i, 1:n)
    end do
  end subroutine print_matrix

  subroutine read_matrix(id, n, matrix)
    integer :: id, i, n
    real, dimension(:,:) :: matrix

    do i = 1, n
       read(id, *) matrix(i, 1:n)
    end do
  end subroutine read_matrix

  subroutine find_solution(A, B, method, n, x)
    integer :: k, i, j, n
    real :: eps, tmp_el
    character :: method
    real, dimension(:,:) :: A
    real, allocatable, dimension(:,:) :: alpha
    ! - матрица альфа -матрица А к которой слева приписали матрицу В. над ней мы и будем проводить преобразования
    real, dimension(:) :: B, X

    allocate(alpha(n,n+1))
    do i = 1,n
       do j = 1,n
          alpha(i,j) = A(i,j)
       end do
       alpha(i,n+1) = B(i)
    end do

    eps = 0.000001
    do i = 1, n
       write(*, *) alpha(i, 1:n+1)
    end do
    select case( method )
       !----------------------------ПРЯМОЙ ХОД ГАУССА----------------------------

    case('g')
       do k = 1,n
          if(abs(alpha(k,k)) .lt. eps) then
             write(*,*) "Attention, value is too close to zero."
          end if

          tmp_el = alpha(k,k)
          forall(j=k:n+1) alpha(k,j) = alpha(k,j)/alpha(k,k)
          forall(i = k+1:n, j = k:n+1) alpha(i,j) = alpha(i,j) - alpha(k,j)*alpha(i,k)
       end do

       do i = 1, n
          write(*, *) alpha(i, 1:n+1)
       end do
       !Вычисляем вектор значений Х
       do i = n, 1, -1
          X(i) = alpha(i,(n+1))
          do j = i+1,n
             X(i) = X(i)- alpha(i,j)*X(j)
          end do
       end do

       write(*,*) X(1:n)
    case default
       write(*,*) "Error, wrong method name.\n"
    end select
  end subroutine find_solution

end module eq_solver
