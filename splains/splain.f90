module splain
  use tri_matr
contains

  subroutine aprox_by_3splains(XYP, res)
    real, dimension(1:100*(size(XYP,dim=1)-1)+1,1:2) :: res
    real, dimension(:,:) :: XYP
    real, dimension(1:size(XYP,dim=1)) :: X, Y, P, C, S, R
    real, dimension(1:size(XYP,dim=1), 1:3) :: A, B, Btran, Q
    real, dimension(1:size(XYP,dim=1), 1:5) :: L, QB
    real :: currx, delta
    integer :: i, n, j


    n = size(XYP,dim=1)

    forall (i=1:n)
       X(i)=XYP(i,1)
       Y(i)=XYP(i,2)
       P(i)=XYP(i,3)
    end forall

    call fill_A_B_Q(X,P,A,B,Q)

    call prepare_eq_set(A,B,Q,Y,L,C)

    call pentadiagonal_matrix_method(L,S,C)

    ! Считаем R - вектор результатов !
    R = 0; QB = 0;
    call tri_transpose(B, Btran)
    call tri_matrix_mult(n, Q, Btran, QB)

    ! Можно тут искать еще !
    R(1) = Y(1) - S(1)*QB(1,3) - S(2)*QB(1,4) - S(3)*QB(1,5)
    R(2) = Y(2) - S(1)*QB(2,2) - S(2)*QB(2,3) - S(3)*QB(2,4) - S(4)*QB(2,5)
    forall(i = 3:n-2) R(i) = Y(i) - S(i-2)*QB(i,1) - S(i-1)*QB(i,2) - S(i)*QB(i,3) - S(i+1)*QB(i,4) - S(i+2)*QB(i,5)
    R(n-1) = Y(n-1) - S(n-1)*QB(n-1,4) - S(n-2)*QB(n-1,3) - S(n-3)*QB(n-1,2) - S(n-4)*QB(n-1,1)
    R(n) = Y(n) - S(n-2)*QB(n,1) - S(n-1)*QB(n,2) - S(n)*QB(n,3)


    ! Тут мы вычисляем значения аппроксимирующей функции на равномерной сетке !
    !                    с числом узлов равным 100n+1                        !

    currx = X(1)
    delta = (X(n) - X(1))/100/(n-1)
    i = 1

    h = X(2) - X(1)
    do j = 1, 100*(n-1)
       if((currx > X(i+1)) .and. (i .lt. n)) then
          i = i+1
          h = X(i+1) - X(i)

          !print *, i, ' ', j
       end if

       t = (currx - X(i))/h
       res(j,1) = currx
       res(j,2) = (1-t)*R(i) + R(i+1)*t - h**2*t*(1-t)/6*((2-t)*S(i) + (1+t)*S(i+1))
       currx = currx + delta
    end do

    res(100*(n-1)+1,1) = currx+delta
    res(100*(n-1)+1,2) = R(i+1)

  end subroutine aprox_by_3splains

  subroutine fill_A_B_Q(X, P, A, B, Q)                                             ! может быть тут
    integer :: i, n
    real, dimension(:) :: X, P
    real, dimension(:,:) :: A,B,Q
    A = 0
    B = 0
    Q = 0
    n = size(X) - 1

    ! Заполняем трехдиагональную матрицу А !

    A(1,2) = 2*(X(2)-X(1)); A(n+1,2) = 2*(X(n+1)-X(n))
    A(2,2) = 2*(X(3)-X(1)); A(2,3) = X(3) - X(2)
    A(n,1) = X(n) - X(n-1); A(n,2) = 2*(X(n+1)-X(n-1));
    A(n+1,2) = 2*(X(n+1) - X(n))
    forall(i = 3:n-1)
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

    forall(i=1:n+1) Q(i,2) = 1/P(i)

  end subroutine fill_A_B_Q

  ! Подсчет матрцы коэффициентов и вектора свободных коэффициентов !
  !      для системы уравнений, решаемой по методу подгонки.       !
  !    В матричном виде система пусть имеет вид LS = C, где        !
  !            L = A + 6*B*Q*Btran, a C = 6*B*Y                    !

  subroutine prepare_eq_set(A,B,Q,Y,L,C)                                     ! В этой процедуре потенциально может быть ошибка
    integer :: i, n
    real, dimension(1:size(Y)) :: Y, BY, C
    real, dimension(1:size(Y), 1:3) :: A, B, Q, Btran, BQ1
    real, dimension(1:size(Y),1:5) :: BQ,  BQBtran, L

    n = size(Y)

    call tri_matrix_mult(n, B, Q, BQ)

    call tri_transpose(B, Btran)

    forall (i=1:n) BQ1(i, 1:3) = BQ(i, 2:4)
    call tri_matrix_mult(n, BQ1, Btran, BQBtran)

    BQBtran = 6*BQBtran
    BQBtran(:,2:4) = BQBtran(:,2:4) + A(:, 1:3)

    L = BQBtran

    BY(1) = B(1,2)*Y(1) + B(1,3)*Y(2)
    do i = 2,n-1
       BY(i) = B(i,1)*Y(i-1) + B(i,2)*Y(i) + B(i,3)*Y(i+1)
    end do
    BY(n) = B(n,2)*Y(n-1) + B(n,3)*Y(n)
    C = BY*6
  end subroutine prepare_eq_set

  subroutine pentadiagonal_matrix_method(L,S,D)
    real, dimension(:) :: D, S
    real, dimension(1:size(D), 1:5) :: L
    real, dimension(-1:size(S)) :: b,a,c,p,q,r, alpha, beta
    integer :: i,n

    n = size(D)

    a(1:n) = L(1:n,3)
    b(1:n-1) = L(2:n,2)
    c(1:n-2) = L(3:n,1)

    a(-1) = 0; a(0) = 0; b(-1) = 0; b(0) = 0; b(n) = 0; c(-1) = 0; c(0) = 0; c(n) = 0;
    c(n-1) = 0; p(-1) = 0; p(0) = 0; q(-1) = 0; q(0) = 0; r(-1) = 0; r(0) = 0;
    alpha(-1) = 0; alpha(0) = 0; beta(-1) = 0; beta(0) = 0;

    do i = 1,n
       beta(i) = b(i-1) - p(i-2)*c(i-2)
       alpha(i) = a(i) - p(i-1)*beta(i) - q(i-2)*c(i-2)
       p(i) = (b(i) - q(i-1)*beta(i))/alpha(i)
       q(i) = c(i)/alpha(i)
       r(i) = (D(i) - r(i-1)*beta(i) - r(i-2)*c(i-2))/alpha(i)
    end do

    ! Ищем наконец вектор решений S !
    S(n) = r(n)
    S(n-1) = r(n-1) - p(n-1)*S(n)

    do i = n-2, 1, -1
       S(i) = r(i) - p(i)*S(i+1) - q(i)*S(i+2)
    end do
  end subroutine pentadiagonal_matrix_method

end module splain
