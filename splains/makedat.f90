module makedat
contains

  !    Генерирует данные для проверки      !
  !         f(x) = x**2 -x + 3             !
  ! Вес допустим равным 3000 в любой точке !

  subroutine make_data_for_5_task(n, x1, x2)
    integer :: i, n

    real :: left, x1, x2, delta

    delta = (x2-x1)/n
    left = x1
    open(10, file = 'data.dat')
    write(10,"('# ',I6)") n

    do i=0,n
       write(10,*) left, sin(left)**2 + 3*log(left), 3000
       left = left + delta
    end do
    close(10)
  end subroutine make_data_for_5_task

end module makedat
