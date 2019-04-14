module eq_solver
contains
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
end module eq_solver