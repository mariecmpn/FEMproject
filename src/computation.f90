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


    subroutine gradient_conjugue(A, L, conv, n)
        integer, intent(in)  :: n ! Dimension des vecteurs et matrices
        real(rp), dimension(n, n), intent(in) :: A ! Matrice du système lineaire
        real(rp), dimension(n), intent(inout) :: L ! en entree: vecteur L second membre du systeme. En sortie: solution du systeme U
        real(rp), dimension(n) :: X, prod ! vecteurs pour la resolution
        real(rp), dimension(n) :: d, gradJ0, gradJ1  ! idem
        real(rp) :: r, rho, conv ! Variables pour le résidu et le pas
        integer :: iter = 0 ! Variable pour le compteur d'itérations
      
        ! Initialisation
        X = 1.E-2
        gradJ0 = matmul(A, X) - L
        d = -gradJ0
        r = 1.
        iter = 0
      
        ! Iterations
        do while ((r > conv) .AND. (iter < 1000))
            iter = iter + 1
            ! Calcul du pas
            prod = matmul(A, d)
            rho = -prod_scal(gradJ0, d, n) / prod_scal(d, prod, n)
            ! Mises a jour
            X = X + rho * d ! solution
            gradJ1 = gradJ0 + rho * matmul(A, d) ! gradient
            ! Mise à jour de la direction de recherche
            d = -gradJ1 + dnrm2(n,gradJ1,1)**2 / dnrm2(n,gradJ0,1)**2 * d ! direction de descente
            r = norme_inf(gradJ1, n) ! norme infinie du gradient
            gradJ0 = gradJ1 
        end do
        L = X ! on enregistre la solution dans L
        write(6,*) 'Nombre d iterations gradient conjugue: ', iter
    end subroutine gradient_conjugue

 
    subroutine prod_matvect_creux(mat, n, X, Y)
        ! Subroutine qui calcule le produit de d'une matrice en stockage creux mat par un vecteur X
        type(mat_creuse), intent(in) :: mat
        integer, intent(in) :: n
        real(rp), dimension(n), intent(in) :: X
        real(rp), dimension(n), intent(out) :: Y
        integer :: i, j, jj, non0
        integer :: c = 0 ! Indicateur de position dans A
        
        Y = 0._rp
        !curseur = 0
        !i = 1
        j = 0
        c = 0
        do i = 1,n
            !non0 = mat%Jposi(i+1) - mat%Jposi(i)
            do j = mat%Jposi(i), mat%Jposi(i+1)-1
                !jj = mat%JvCell(j)
                Y(i) = Y(i) + mat%Tmat(j)*X(mat%JvCell(j))
            end do
        end do
        !DO WHILE (curseur < mat%NCoefMat)
        !    non0 = mat%Jposi(i+1) - mat%Jposi(i)
        !    do ii = 1,non0
        !        curseur = curseur + 1
        !        j = mat%JvCell(curseur) ! Colonne du terme non-nul où est le curseur
        !        Y(i) = Y(i) + mat%Tmat(curseur) * X(j)
        !    end do
        !    i = i+1
        !END DO
        
    end subroutine prod_matvect_creux

  

    SUBROUTINE Gradient_conjugue_CREUX (A, Y, conv, dim_mat)
        IMPLICIT NONE
        INTEGER, INTENT(in)             ::  dim_mat
        type(mat_creuse), INTENT(in)         ::  A
        REAL(rp), DIMENSION(dim_mat), INTENT(inout)  ::  Y
        REAL(rp), DIMENSION(dim_mat) ::  X
        REAL(rp), DIMENSION(dim_mat) :: d, gradJ0, gradJ1, prod
        REAL(rp) :: rho, r
        INTEGER :: iter
        real(rp), intent(in) :: conv
        
        X = 0.
        CALL prod_matvect_creux(A, dim_mat, X, prod) 
        gradJ0 = prod - Y
        d = -gradJ0
        r = 1.
        iter = 0
        DO WHILE ((r > conv) .AND. (iter < 1000))
            iter = iter + 1
            CALL prod_matvect_creux(A, dim_mat, d, prod) ! tmp = A.d
            rho = -prod_scal(gradJ0, d, dim_mat) / prod_scal(d, prod, dim_mat)
            X = X + rho * d
            gradJ1 = gradJ0 + rho * prod
            d = -gradJ1 + sum(gradJ1**2) / sum(gradJ0**2) * d
            r = norme_inf(gradJ1, dim_mat)
            gradJ0 = gradJ1
        END DO
        Y = X
        write(6,*) 'Nombre d iterations gradient conjugue: ', iter
    END SUBROUTINE Gradient_Conjugue_CREUX

end module computation