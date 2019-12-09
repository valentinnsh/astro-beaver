module file_ops
  implicit none
contains

  function get_file_name(n) result(filename)
    integer(4), intent(in) :: n
    character(len=14) :: filename
    character(len=12) :: format_string

    if (n < 10) then
       format_string = "(A4, I1, A4)"
    else
       format_string = "(A4, I2, A4)"
    endif
    write (filename,format_string) "quad", n, ".dat"

  end function get_file_name

  recursive subroutine quicksort(a, first, last)
    implicit none
    real(8)  a(*), x, t
    integer(4) :: i, j, first, last

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do

       do while (x < a(j))
          j=j-1
       end do

       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do

    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort

  function f(x)
    real(8) :: x, f
    f = x + 1
  end function f

end module file_ops
