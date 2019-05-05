module splains
  use tri_matr
contains

  subroutine aprox_by_3splains(XYP, res)
    real, dimension(:,1:2) :: res
    real, dimension(:,:) :: XYP
    real, dimension(1:size(XYP,dim=2)) :: X, Y, P, C, S, R
    real, dimension(1:size(XYP,dim=2), 1:3) :: A, B, Btran, Q
    real, dimension(1:size(XYP,dim=2), 1:5) :: L, QB
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


end module splains
