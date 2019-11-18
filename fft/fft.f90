module fft
  implicit none
contains
  
  recursive subroutine bit_reversal_sorting(X, Y, Wn) ! Реверсивная сортировка
    complex, dimension(0:), intent(in) :: X
    complex, dimension(0:size(X) - 1), intent(out) :: Y
    complex, dimension(0:), intent(in) :: Wn
    integer :: i, N

    N = size(X)

    if (N>2) then
       call bit_reversal_sorting(X(0:N/2-1) + X(N/2:N-1),               Y(0:N-2:2), Wn(0:N-1:2))
       call bit_reversal_sorting(Wn(0:N/2-1)*(X(0:N/2-1) - X(N/2:N-1)), Y(1:N-1:2), Wn(0:N-1:2))
    else
       Y(0) = X(0) + X(1)
       Y(1) = X(0) - X(1)
    endif
    
  end subroutine bit_reversal_sorting

  
  subroutine W_negative(W, N)
    integer :: q, N
    complex, dimension(0:) :: W

    forall (q=0:N-1) W(q) = exp(-2.0*(4*atan(1.0))*cmplx(0,1)*q/N)

  end subroutine  W_negative

  subroutine W_positive(W, N)
    integer :: q, N
    complex, dimension(0:) :: W

    forall (q=0:N-1) W(q) = exp(2.0*(4*atan(1.0))*cmplx(0,1)*q/N)

  end subroutine  W_positive

  
end module fft
