program tester
  use makedat
  implicit none

  integer :: i, j, n
  real, allocatable, dimension(:,:) :: XYP

  call make_data_for_5_task(10,0.0,100.0)

  open(10, file = 'data.dat')
  read(10,'(2x,I6)') n
  allocate(XYP(1:n+1,3))

  do i = 1, n+1
     read(10,*) XYP(i,1:3)
  end do


  close(10)
  deallocate(XYP)
end program tester
