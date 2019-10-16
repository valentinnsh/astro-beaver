program tester
  use makedat
  use splain
  implicit none

  integer :: i, n
  real, allocatable, dimension(:,:) :: XYP
  real, allocatable, dimension(:,:) :: res



  call make_data_for_5_task(5,0.0,10.0)

  open(10, file = 'data.dat')
  open(20, file = 'result.dat')
  read(10,'(2x,I6)') n
  allocate(XYP(n+1,3))
  allocate(res(1:100*n+1,2))

  do i = 1, n+1
     read(10,*) XYP(i,1:3)
  end do

  res = 0

  call aprox_by_3splains(XYP,res)


  do i = 1,100*n+1
     write(20,*) res(i,1), res(i,2)
  end do

  close(10)
  close(20)
  deallocate(XYP)
  deallocate(res)
end program tester
