module tri_makedata
contains
  subroutine make_data(n)
    integer i, j, n
    real tmp
    real, allocatable, dimension(:,:) :: matrix
    open(10, file = 'data1.dat')
    write(10,"('# ',I6)") n
    allocate(matrix(n,3))
    call RANDOM_NUMBER(matrix)

    write(10,*) matrix(1, 1:2)
    i = 2
    do while (i .le. n-1)
       write(10,*) matrix(i, 1:3)
       i = i+1
    end do
    write(10,*) matrix(n, 1:2)

    close(10)
    deallocate(matrix)

    open(10, file = 'data2.dat')
    write(10,"('# ',I6)") n
    allocate(matrix(n,3))
    call RANDOM_NUMBER(matrix)

    write(10,*) matrix(1, 1:2)
    i = 2
    do while (i .le. n-1)
       write(10,*) matrix(i, 1:3)
       i = i +1
    end do
    write(10,*) matrix(n, 1:2)
    
    close(10)
    deallocate(matrix)
    
  end subroutine make_data
end module tri_makedata
