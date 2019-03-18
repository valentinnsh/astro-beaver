module matrix
contains
  subroutine native_mult(n, matrix_A, matrix_B, matrix_C) ! "Лобовой" алгоритм перемножения матриц
    integer i, j, k, n
    real,  dimension(:, :) :: matrix_A, matrix_B, matrix_C
    real  start, finish
    call cpu_time(start)

    do i = 1,n
       do j = 1,n
          do k = 1,n
             matrix_C(i,j) = matrix_C(i,j) + matrix_A(i,k)*matrix_B(k,j)
          end do
       end do
    end do
    

    call cpu_time(finish)
    write(*,*) finish - start
  end subroutine native_mult

  subroutine print_matrix(id, n, matrix)
    integer :: id, i, n
    real, dimension(:,:) :: matrix
    do i = 1, n
       write(id, *) matrix(1:n, i)
    end do
  end subroutine print_matrix

  subroutine read_matrix(id, n, matrix)
    integer :: id, i, n
    real, dimension(:,:) :: matrix

    do i = 1, n
       read(id, *) matrix(1:n, i)
    end do
  end subroutine read_matrix
  
end module matrix
