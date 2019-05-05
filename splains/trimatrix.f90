module tri_matr
contains

  subroutine tri_transpose(A, tranA)
    real, dimension(:,:) :: A, tranA
    integer ::i,n

    tranA = 0
    tranA(0,2) = A(0,2); tranA(0,3) = A(1,1)

    forall(i = 1:n-1)
       tranA(i,1) = A(i-1,3)
       tranA(i,2) = A(i,2)
       tranA(i,3) = A(i+1,1)
    end forall

    tranA(n,1) = A(n-1,3); tranA(n,2) = A(n,2)
  end subroutine tri_transpose

  subroutine tri_matrix_mult(n, matrix_A, matrix_B, matrix_C)
    integer i, j, k, n
    real,  dimension(:, :) :: matrix_A, matrix_B
    real, dimension(:,:) :: matrix_C
    real :: start, finish

    ! Считаем первый столбик
    do i=3,n
       matrix_C(i,1)=matrix_A(i,1)*matrix_B(i-1,1)
    enddo

    !Второй столбец
    do i=2,n
       matrix_C(i,2)=matrix_A(i,1)*matrix_B(i-1,2)+matrix_A(i,2)*matrix_B(i,1)
    enddo

    !Третий столбик
    matrix_C(1,3)=matrix_A(1,2)*matrix_B(1,2)+matrix_A(1,3)*matrix_B(2,1)
    do i=2,n-1
       matrix_C(i,3)=matrix_A(i,1)*matrix_B(i-1,3)+matrix_A(i,2)*matrix_B(i,2)+matrix_A(i,3)*matrix_B(i+1,1)
    enddo

    matrix_C(n,3)=matrix_A(n,1)*matrix_B(n-1,3)+matrix_A(n,2)*matrix_B(n,2)

    !Четвертый столбец
    do i=1,n-2
       matrix_C(i,4)=matrix_A(i,2)*matrix_B(i,3)+matrix_A(i,3)*matrix_B(i+1,2)
    enddo

    matrix_C(n-1,4)=matrix_A(n-1,2)*matrix_B(n-1,3)+matrix_A(n-1,3)*matrix_B(n,2)

    !Пятый столбец
    do i=1,n-2
       matrix_C(i,5)=matrix_A(i,3)*matrix_B(i+1,3)
    enddo

    write(*,*) matrix_C(1, 3)
  end subroutine tri_matrix_mult

end module tri_matr
