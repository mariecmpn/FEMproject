module computation
    use numerics
    IMPLICIT NONE

    !--------------------------------------------!
    !               DANS CE MODULE:
    ! - Fonctions et subroutines utiles pour la
    ! resolution du systeme lineaire
    !--------------------------------------------!

    real(rp), external :: dnrm2 ! fonction de BLAS pour la norme L2

    contains

    real(rp) function norme_inf(X, n)
        ! fonction qui calcule la norme infinie d'un vecteur
        integer :: n
        real(rp), dimension(n) :: X
        integer :: i
        real(rp) :: m
        m = 0._rp
        do i = 1,n
            m = max(m, abs(X(i)))
        end do
        norme_inf = m
    end function norme_inf

    real(rp) function dot_trans(X, Y, dim_mat)
        ! fonction qui calcule le produit X^T * Y, ou X et Y sont des veteurs de dimension dim_mat
        integer :: dim_mat ! dimension des vecteurs
        real(rp), dimension(dim_mat) :: X,Y ! vecteurs dont on souhaite connaitre le produit
        integer :: i
        real(rp) :: d
        d = 0._rp
        do i = 1,dim_mat
            d = d+X(i)*Y(i)
        end do
        dot_trans = d
    end function dot_trans

    function dot_mat_vec(Mat, Vec, dim_mat)
        integer :: dim_mat
        real(rp), dimension(dim_mat) :: dot_mat_vec
        real(rp), dimension(dim_mat, dim_mat) :: Mat
        real(rp), dimension(dim_mat) :: Vec
        integer :: i,j
        dot_mat_vec(:) = 0._rp
        do j = 1,dim_mat
            do i = 1,dim_mat
            dot_mat_vec(i) = dot_mat_vec(i) + Mat(i,j)*Vec(j)
            end do
        end do
    end function dot_mat_vec

    real(rp) function prod_scal(vec1, vec2, dim_mat)
        integer :: dim_mat
        real(rp), dimension(dim_mat) :: vec1, vec2
        integer :: i
        prod_scal = 0._rp
        do i = 1,dim_mat
            prod_scal = prod_scal+ vec1(i)*vec2(i)
        end do
    end function prod_scal


    !subroutine gradient_conjugue(A, L, conv, dim_mat)
        ! subroutine pour la resolution d'un systeme lineaire Ax = L par une algorithme de gradient conjugue
    !    integer, intent(in) :: dim_mat ! dimension du vecteur et de la matrice
    !    real(rp), dimension(dim_mat,dim_mat), intent(in) :: A ! matrice A
    !    real(rp), dimension(dim_mat), intent(inout) :: L ! en entree: vecteur L second membre. En sortie: solution du systeme x
    !    real(rp), intent(in) :: conv ! precision de convergence souhaitee
        !real(rp), dimension(dim_mat) :: r_k, x_k, r_k1, p_k ! vecteurs residus r_k et r_k1, vecteur direction de descente p_k, vecteur itere x_k
    !    real(rp), dimension(dim_mat) :: Apk ! vecteur pour les multiplications entre A et un vecteur
    !    real(rp), dimension(dim_mat) :: grad, grad1, d_k, x_k, x_k1
    !    integer :: nb_iter, i ! nombre d'iterations
    !    real(rp) :: rho_k, beta, p ! coef alpha_k et beta_k dans la boucle while 

        ! initialisation
    !    x_k(:) = 1.D-2 ! on choisit x_0 de facon arbitraire
        !x_k(:) = 0._rp
    !    Apk = matmul(A,x_k)
        !Apk = dot_mat_vec(A,x_k,dim_mat) ! A*x_0
        !r_k(:) = Apk(:)-L(:) ! r_0 = L-A*x_0
        !p_k(:) = r_k(:)
        !write(6,*) (r_k(i), i = 1,dim_mat)
    !    grad(:) = Apk(:) - L(:)
    !    d_k(:) = -grad(:)
    !    nb_iter = 0

        ! iterations
    !    do while (nb_iter<=1000 .AND. norme_inf(grad, dim_mat)>conv) ! tant qu'on n'a pas depasse le nombre d'iterations max
            !Apk = matmul(A, p_k) 
    !        Apk = matmul(A,d_k) ! A*p_k
    !        write(6,*) (Apk(i), i = 1,dim_mat)
    !        p = prod_scal(Apk, d_k, dim_mat)
    !        write(6,*) p
    !        rho_k = - prod_scal(grad, d_k, dim_mat) / p
    !        x_k1(:) = x_k(:)
    !        x_k(:) = x_k1(:) + rho_k*d_k(:)
    !        grad1(:) = grad(:)
    !        grad(:) = grad1(:) + rho_k*dot_mat_vec(A,d_k,dim_mat)
    !       beta = dnrm2(dim_mat,grad,1)**2 / dnrm2(dim_mat,grad1,1)**2
    !        d_k(:) = grad(:) + beta*d_k(:)
    !        nb_iter = nb_iter+1
            !alpha = dot_trans(r_k, r_k, dim_mat) / dot_trans(p_k, Apk, dim_mat) ! r_k^T*r_k / p_k^T*A*p_k
            !x_k(:) = x_k(:) + alpha*p_k(:) 
            !r_k1(:) = r_k(:) ! On garde en memoire r_k dans r_k1 si besoin de calculer beta
            !r_k(:) = r_k(:) - alpha*Apk(:)
            !if (norme_inf(r_k, dim_mat)<=conv) then ! si r_k est assez petit, on sort de la boucle 
            !    exit
            !else ! sinon on initialise beta et p_k pour la prochaine iteration
                !beta = dot_trans(r_k, r_k, dim_mat) / dot_trans(r_k1, r_k1, dim_mat)
                !p_k(:) = r_k(:) + beta*p_k(:)
                !nb_iter = nb_iter+1
            !end if
    !    end do
    !    write(6,*) 'Nombre d iterations gradient conjugue: ', nb_iter
    !    L(:) = x_k(:) ! la solution du syteme est x_k+1
    !end subroutine gradient_conjugue


    SUBROUTINE Gradient_conjugue(A, L, conv, n)

        IMPLICIT NONE
      
        INTEGER, INTENT(IN)  :: n             ! Dimension des vecteurs et matrices
        REAL(rp), DIMENSION(n, n), INTENT(IN)   :: A  ! Matrice du système linéaire
        REAL(rp), DIMENSION(n), INTENT(INOUT)  :: L  ! Vecteur solution
        REAL(rp), DIMENSION(n) :: X, prod
        REAL(rp), DIMENSION(n) :: d, gradJ0, gradJ1  ! Vecteurs auxiliaires
        REAL(rp) :: r, rho, conv          ! Variables pour le résidu et le pas
        INTEGER :: iter = 0                        ! Variable pour le compteur d'itérations
      
        ! Initialisation du vecteur solution X à zéro
        X = 1.E-2
      
        ! Calcul du gradient initial de J(X):=1/2(A.X).X - Y.X
        gradJ0 = matmul(A, X) - L
        d = -gradJ0
        r = 1.
        iter = 0
      
        ! Boucle principale de la méthode du gradient conjugué
        DO WHILE ((r > conv) .AND. (iter < 1000))
          iter = iter + 1
          ! Calcul du pas
          prod = matmul(A, d)
          rho = -prod_scal(gradJ0, d, n) / prod_scal(d, prod, n)
          ! Mise à jour du vecteur solution X
          X = X + rho * d
          ! Mise à jour du gradient
          gradJ1 = gradJ0 + rho * matmul(A, d)
          ! Mise à jour de la direction de recherche
          d = -gradJ1 + dnrm2(n,gradJ1,1)**2 / dnrm2(n,gradJ0,1)**2 * d
          ! Mise à jour du résidu
          r = maxval(gradJ1)
          gradJ0 = gradJ1
        END DO
        L = X
        write(6,*) 'Nombre d iterations gradient conjugue: ', iter
    END SUBROUTINE Gradient_conjugue


  ! routine pour calculer le produit d une matrice creuse A par un vecteur  
  SUBROUTINE PRODUIT_MxV_CREUX (Mat, n, X, Y)
    IMPLICIT NONE ! Mat . X = Y
    type(mat_creuse), INTENT(in) :: Mat
    INTEGER, INTENT(in) :: n
    REAL(rp), DIMENSION(n), INTENT(in) :: X
    REAL(rp), DIMENSION(n), INTENT(out) :: Y
    INTEGER :: i = 0, j = 0, ii, non0
    INTEGER :: curseur = 0 ! Indicateur de position dans A
    
    Y = 0.
    curseur = 0
    i = 1
    DO WHILE (curseur < Mat%NCoefMat)
        non0 = Mat%Jposi(i+1) - Mat%Jposi(i)
        do ii = 1,non0
            curseur = curseur + 1
            j = Mat%JvCell(curseur) ! Colonne du terme non-nul où est le curseur
            Y(i) = Y(i) + Mat%Tmat(curseur) * X(j)
        end do
        i = i+1
    END DO
    
  END SUBROUTINE PRODUIT_MxV_CREUX

  

  SUBROUTINE Gradient_conjugue_CREUX (A, Y, conv, n)
    IMPLICIT NONE
    INTEGER, INTENT(in)             ::  n
    type(mat_creuse), INTENT(in)         ::  A
    REAL(rp), DIMENSION(n), INTENT(inout)  ::  Y
    REAL(rp), DIMENSION(n) ::  X
    REAL(rp), DIMENSION(n) :: d, gradJ0, gradJ1, tmp
    REAL(rp) :: rho = 0., residu = 0.
    INTEGER :: iteration = 0
    real(rp), intent(in) :: conv
    
    X = 0.
    CALL PRODUIT_MxV_CREUX(A, n, X, tmp) ! tmp = A.X
    gradJ0 = tmp - Y
    d = -gradJ0
    residu = 1.
    iteration = 0
    DO WHILE ((residu > conv) .AND. (iteration < 1000))
       iteration = iteration + 1
       CALL PRODUIT_MxV_CREUX(A, n, d, tmp) ! tmp = A.d
       rho = -dot_product(gradJ0, d) / dot_product(d, tmp)
       X = X + rho * d
       gradJ1 = gradJ0 + rho * tmp
       d = -gradJ1 + sum(gradJ1**2) / sum(gradJ0**2) * d
       residu = maxval(gradJ1)
       gradJ0 = gradJ1
    END DO
    Y = X
  END SUBROUTINE Gradient_Conjugue_CREUX

end module computation