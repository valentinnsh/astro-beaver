module matrixop
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
end module matrixop
