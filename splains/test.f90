program tester
  use makedat
  implicit none

  integer :: i, j, n

  call make_data_for_5_task(10,0.0,100.0)
  open(10, file = 'data.dat')

  close(10)
end program tester
