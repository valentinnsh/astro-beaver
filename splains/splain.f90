module splains
  use tri_matr
contains

  subroutine aprox_by_3splains(XYP, res)
    real, dimension(:,1:2) :: res
    real, dimension(:,:) :: XYP
    real, dimension(1:size(XYP,dim=2)) :: X, Y, P, C, S, R
    real, dimension(1:size(XYP,dim=2), 1:3) :: A, B, Btran, Q
    real, dimension(1:size(XYP,dim=2), 1:5) :: L, QB
    real :: currx, delta
    integer :: i, n


    n = size(XYP,dim=2) - 1

    forall (i=1:n+1)
       X(i)=XYP(i,1)
       Y(i)=XYP(i,2)
       P(i)=XYP(i,3)
    end forall

    call fill_A_B_Q(X,Y,P,A,B,Q)
    call prepare_eq_set(A,B,Q,Y,L,C)
    call pentadiagonal_matrix_method(L,S,C)

    ! Считаем R - вектор результатов !
    R = 0; QB = 0;
    call tri_transpose(B, Btran)
    call tri_matrix_mult(n+1, Q, Btran, QB)

    R(1) = Y(1) - S(1)*QB(1,3) - S(2)*QB(1,4) - S(3)*QB(1,5)
    R(2) = Y(2) - S(1)*QB(2,2) - S(2)*QB*(2,3) - S(3)*QB(2,4) - S(4)*QB(2,5)
    forall(i = 3:n-1) R(i) = Y(i) - S(i-2)*QB(i,1) - S(i-1)*QB(i,2) - S(i)*QB(i,3) - S(i+1)*QB(i,4) - S(i+2)*QB(i,5)
    R(n) = Y(n) - S(n+1)*QB(n,4) - S(n)*QB(n,3) - S(n-1)*QB(n,2) - S(n-2)*QB(n,1)
    R(n+1) = Y(n+1) - S(n-1)*QB(n+1,1) - S(n)*QB(n+1,2) - S(n+1)*QB(n+1,3)

    ! Тут мы вычисляем значения апроксимирующей функции на равномерной сетке !
    !                    с числом узлов равным 100n+1                        !

    currx = X(1)
    delta = (X(n+1) - X(1))/100/n
    do j = 1, 100*n+1
       if(currx > X(i)) i = i+1
       h = X(i+1) - X(i)
       t = (currx - X(i))/h
       res(j,1) = currx
       res(j,2) = (1-t)*R(i) + R(i+1)*t - h**2*t*(1-t)/6*((2-t)*S(i) + (1+t)*S(i+1))
       currx = currx + delta
    end do

  end subroutine aprox_by_3splains

  subroutine fill_A_B_Q(X,Y,P,A,B,Q)
    integer :: i, n
    real, dimension(:) :: X, Y, P
    real, dimension(:,:) :: A,B,Q
    A = 0
    B = 0
    Q = 0
    n = size(X) - 1

    ! Заполняем трехдиагональную матрицу А !

    A(1,2) = 2*(X(2)-X(1)); A(n+1,2) = 2*(X(n+1)-X(n))
    A(2,2) = 2*(X(3)-X(1)); A(2,3) = X(3) - X(2)
    A(n,1) = X(n) - X(n-1); A(n,2) = 2*(X(n+1)-X(n-1))
    forall(i = 3:n)
       A(i,1) = X(i) - X(i-1)
       A(i,2) = 2*(X(i+1) - X(i-1))
       A(i,3) = X(i+1) - X(i)
    end forall

    ! Ззаполняем трехдиагональную матрицу В !

    forall(i = 2:n)
       B(i,1) = 1/(X(i) - X(i-1))
       B(i,2) = -1/(X(i) - X(i-1)) - 1/(X(i+1) - X(i))
       B(i,3) = 1/(X(i+1) - X(i))
    end forall

    ! Заполняем трехдиагональную матрицу Q !

    forall(i=1,n+1) Q(i,2) = 1/P(i)

  end subroutine fill_A_B_Q

  ! Подсчет матрцы коэффициентов и вектора свободных коэффициентов !
  !      для системы уравнений, решаемой по методу подгонки.       !
  !    В матричном виде система пусть имеет вид LS = C, где        !
  !            L = A + 6*B*Q*Btran, a C = 6*B*Y                    !

  subroutine prepare_eq_set(A,B,Q,Y,L,C)
    integer :: i, n
    real, dimension(:) :: Y, P, C
    real, dimension(:,:) :: A, B, Q
    real, dimension(1:size(Y),1:5) :: BQ, BQBtran

    n = size(Y) - 1
    call tri_matrix_mult(n-1, B, Q, BQ)

  end subroutine prepare_eq_set

  subroutine pentadiagonal_matrix_method(L,S,C)
    real, dimension(:) :: C, S
    real, dimension(:, 1:5) :: L
    real, dimension(1:size(S)) :: b,a,p,q,r
    integer :: i,n

    n = size(C) - 1

    b(1) = 0
    a(1) = L(1,3)
    p(1) = L(2,2)/a(1)
    q(1) = C(1)
    r(1) = C(1)/a(1)

    b(2) = L(2,2)
    a(2) = L(2,3) - p(1)*b(2)
    p(2) = (L(3,2) - q(1)*b(2))/a(2)
    q(2) = L(4,1)/a(2)
    r(2) = (C(2)-r(1)*b(2))/a(2)

    do i = 3,n-1
       b(i) = L(i,2) - p(i-2)*L(i,1)
       a(i) = L(i,3) - p(i-1)*b(i) - q(i-2)*L(i,1)
       p(i) = (L(i+1,2) - q(i-1)*b(i))/a(i)
       q(i) = L(i+2,1)/a(i)
       r(i) = (C(i) - r(i-1)*b(i) - r(i-2)*L(i,1))/a(i)
    end do

    b(n) = L(n,2) - p(n-2)*L(n,1)
    a(n) = L(n,3) - p(n-1)*b(n) - q(n-1)*L(n,1)
    p(n) = (L(n+1,2) - q(n-1)*b(n))/a(n)
    r(n) = (C(n) - r(n-1)*b(n) - r(n-2)*L(n,1))/a(n)
    q(n) = 0

    b(n+1) = L(n+1,2) - p(n-1)*L(n+1,1)
    a(n+1) = L(n+1,3) - p(n)*b(n+1) - q(n)*L(n+1,1)
    p(n+1) = (0 - q(n)*b(n+1))/a(n+1)
    q(n+1) = 0
    r(n+1) = (C(n+1) - r(n)*b(n+1) - r(n-1)*L(n+1,1))/a(n+1)

    ! Ищем наконец вектор решений S !
    S(n+1) = r(n+1)
    S(n) = r(n) - p(n)*S(n+1)

    do i = n-1, 1, -1
       S(i) = r(i) - p(i)*S(i+1) - q(i)*S(i+2)
    end do
  end subroutine pentadiagonal_matrix_method

end module splains
