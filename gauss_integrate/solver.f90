module eq_solver
contains
  subroutine print_matrix(id, n, matrix)
    integer(4) :: id, i, n
    real(8), dimension(:,:) :: matrix
    do i = 1, n
       write(id, *) matrix(i, 1:n)
    end do
  end subroutine print_matrix

  subroutine read_matrix(id, n, matrix)
    integer(4) :: id, i, n
    real(8), dimension(:,:) :: matrix

    do i = 1, n
       read(id, *) matrix(i, 1:n)
    end do
  end subroutine read_matrix

  subroutine find_solution(A, B, method, n, X)
    integer(4) :: k, i, j, n
    real(8) :: eps, tmp_el
    character :: method
    real(8), dimension(:,:) :: A
    real(8), dimension(n+1) :: tmp_X
    real(8), allocatable, dimension(:,:) :: alpha
    integer(4), allocatable, dimension(:,:) :: swap_logs
    ! - матрица альфа -матрица А к которой слева приписали матрицу В. над ней мы и будем проводить преобразования
    real(8), dimension(:) :: B, X

    allocate(alpha(n,n+1))
    do i = 1,n
       do j = 1,n
          alpha(i,j) = A(i,j)
       end do
       alpha(i,n+1) = B(i)
    end do

    eps = 0.0000000001_8

    select case( method )
       !----------------------------Гаусс----------------------------
    case('g')
       do k = 1,n
          if(abs(alpha(k,k)) .lt. eps) then
             write(*,*) "Attention, value is too close to zero."
          end if

          tmp_el = alpha(k,k)
          forall(j=k:n+1) alpha(k,j) = alpha(k,j)/tmp_el
          forall(i = k+1:n, j = k:n+1) alpha(i,j) = alpha(i,j) - alpha(k,j)*alpha(i,k)
       end do

       !Вычисляем вектор значений Х
       do i = n, 1, -1
          X(i) = alpha(i,(n+1))
          do j = i+1,n
             X(i) = X(i)- alpha(i,j)*X(j)
          end do
       end do
       !-------Жордан---------------------------------
    case('j')
       do k = 1,n
          if(abs(alpha(k,k)) .lt. eps) then
             write(*,*) "Attention, value is too close to zero."
          end if

          tmp_el = alpha(k,k)
          forall(j=k:n+1) alpha(k,j) = alpha(k,j)/tmp_el
          forall(j = k:n+1, i =1:n, i .ne. k) alpha(i,j) = alpha(i,j) - alpha(k,j)*alpha(i,k)
       end do

       do i = 1,n
          X(i) = alpha(i,n+1)
       end do
       !-----------Выбор ведущего элемента --------------
    case('i')
       allocate(swap_logs(n,2))
       swap_logs = 0
       do k = 1,n-1

          if(abs(maxval(alpha(k:n,k:n))) .gt. abs(minval(alpha(k:n,k:n)))) then
             swap_logs(k,1:2) = k-1 + maxloc(alpha(k:n,k:n))
          else
             swap_logs(k,1:2) = k-1 + minloc(alpha(k:n,k:n))
          end if
          tmp_X(1:n) = alpha(1:n,swap_logs(k,2))
          alpha(1:n, swap_logs(k,2)) = alpha(1:n,k)
          alpha(1:n,k) = tmp_X(1:n)

          tmp_X(1:n+1) = alpha(swap_logs(k,1), 1:n+1)
          alpha(swap_logs(k,1), 1:n+1) = alpha(k, 1:n+1)
          alpha(k, 1:n+1) = tmp_X(1:n+1)

          forall(j=k:n+1) alpha(k,j) = alpha(k,j)/alpha(k,k)
          forall(i = k+1:n, j = k:n+1) alpha(i,j) = alpha(i,j) - alpha(k,j)*alpha(i,k)

       end do

       alpha(n,n:n+1) = alpha(n,n:n+1)/alpha(n,n)

       !Вычисляем вектор значений Х
       do i = n, 1, -1
          X(i) = alpha(i,(n+1))
          do j = i+1,n
             X(i) = X(i)- alpha(i,j)*X(j)
          end do
       end do

       do i = n-1, 1, -1
          tmp_X(i) = X(swap_logs(i,1))
          X(swap_logs(i,1)) = X(i)
          X(i) = tmp_X(i)
       end do
       deallocate(swap_logs)
    case default
       write(*,*) "Error, wrong method name.\n"
    end select
    deallocate(alpha)
  end subroutine find_solution

end module eq_solver
