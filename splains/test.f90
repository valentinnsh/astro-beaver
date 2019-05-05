program tester
  use makedat
  implicit none

  integer :: i, j, n
  real, allocatable, dimension(:,:) :: XYP

  call make_data_for_5_task(10,0.0,100.0)

  open(10, file = 'data.dat')
  read(10,'(2x,I6)') n
  allocate(XYP(0:n,3))

  do i = 0,n
     read(10,*) XYP(i,1:3)
  end do


  close(10)
  deallocate(XYP)
end program tester
