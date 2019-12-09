program tester
  use horner_and_bernoulli
  use legendre
  use file_ops
  use gauss
  implicit none

  character(len=14) :: filename
  character(len=2) :: n1
  real(8), allocatable, dimension(:) :: A, t
  real(8), parameter :: lim_a = 0.0_8, lim_b = 1.0_8
  real(8) :: res
  integer(4) :: i, n

  interface
     function func(x)
       real(8) :: x, func
     end function func
  end interface

  call getarg(1, n1)
  read(n1, *) n
  allocate(A(1:n), t(1:n))
  filename = get_file_name(n)

  !           _,'|             _.-''``-...___..--';)
  !         /_ \'.      __..-' ,      ,--...--'''
  !        <\    .`--'''       `     /'
  !         `-';'               ;   ; ;
  !   __...--''     ___...--_..'  .;.'
  !  (,__....----'''       (,..--''

  t = legendre_roots(n)

  A = gauss_quad_coef(n, t)

  open(10, file=filename)
  do i = 1, n
     write(10, '(e15.7,3x,e15.7)') A(i), t(i)
  enddo
  close(10)

  res = integration(lim_a, lim_b, n, f)

  write(*,*) res

end program tester
