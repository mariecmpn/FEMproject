module computation
    use numerics
    IMPLICIT NONE

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


    subroutine gradient_conjugue(A, L, conv, dim_mat)
        ! subroutine pour la resolution d'un systeme lineaire Ax = L par une algorithme de gradient conjugue
        integer, intent(in) :: dim_mat ! dimension du vecteur et de la matrice
        real(rp), dimension(dim_mat,dim_mat), intent(in) :: A ! matrice A
        real(rp), dimension(dim_mat), intent(inout) :: L ! en entree: vecteur L second membre. En sortie: solution du systeme x
        real(rp), intent(in) :: conv ! precision de convergence souhaitee
        real(rp), dimension(dim_mat) :: r_k, x_k, r_k1, p_k ! vecteurs residus r_k et r_k1, vecteur direction de descente p_k, vecteur itere x_k
        real(rp), dimension(dim_mat) :: Apk ! vecteur pour les multiplications entre A et un vecteur
        integer :: nb_iter,i ! nombre d'iterations
        real(rp) :: alpha, beta ! coef alpha_k et beta_k dans la boucle while 

        ! initialisation
        x_k(:) = 1.D-2 ! on choisit x_0 de facon arbitraire
        !Apk = matmul(A,x_k)
        Apk = dot_mat_vec(A,x_k,dim_mat) ! A*x_0
        write(6,*) 'OK'
        r_k(:) = L(:)-Apk(:) ! r_0 = L-A*x_0
        p_k(:) = r_k(:)
        !write(6,*) (r_k(i), i = 1,dim_mat)
        nb_iter = 0

        ! iterations
        do while (nb_iter<=1000) ! tant qu'on n'a pas depasse le nombre d'iterations max
            !Apk = matmul(A, p_k) 
            Apk = dot_mat_vec(A,p_k,dim_mat) ! A*p_k
            write(6,*) (Apk(i), i = 1,dim_mat)
            write(6,*) dot_trans(p_k, Apk, dim_mat)
            alpha = dot_trans(r_k, r_k, dim_mat) / dot_trans(p_k, Apk, dim_mat) ! r_k^T*r_k / p_k^T*A*p_k
            x_k(:) = x_k(:) + alpha*p_k(:) 
            r_k1(:) = r_k(:) ! On garde en memoire r_k dans r_k1 si besoin de calculer beta
            r_k(:) = r_k(:) - alpha*Apk(:)
            if (norme_inf(r_k, dim_mat)<=conv) then ! si r_k est assez petit, on sort de la boucle 
                exit
                write(6,*) 'OK'
            else ! sinon on initialise beta et p_k pour la prochaine iteration
                beta = dot_trans(r_k, r_k, dim_mat) / dot_trans(r_k1, r_k1, dim_mat)
                p_k(:) = r_k(:) + beta*p_k(:)
                nb_iter = nb_iter+1
            end if
        end do
        write(6,*) 'Nombre d iterations gradient conjugue: ', nb_iter
        L(:) = x_k(:) ! la solution du syteme est x_k+1
    end subroutine gradient_conjugue


end module computation