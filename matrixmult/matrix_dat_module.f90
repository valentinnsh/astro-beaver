module matrix_dat
use matrix
contains
  subroutine make_data(n)
    integer i, j, n
    real tmp
    real, allocatable, dimension(:,:) :: matrix

    open(10, file = 'data1.dat')
    write(10,"('# ',I6)") n
    allocate(matrix(n,n))
    call RANDOM_NUMBER(matrix)
    call print_matrix(10, n, matrix)
    close(10)

    open(10, file = 'data2.dat')
    write(10,"('# ',I6)") n
    call RANDOM_NUMBER(matrix)
    call print_matrix(10, n, matrix)
    close(10)
    
    deallocate(matrix)
  end subroutine make_data
end module matrix_dat
