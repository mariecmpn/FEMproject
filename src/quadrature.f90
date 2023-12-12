module quadrature_P1
    use numerics
    IMPLICIT NONE

    contains

    !--------------------------------------------!
    !             Fonctions exactes
    !--------------------------------------------!

    real(rp) function f(coor)
        real(rp), dimension(2) :: coor
        f = -2._rp * (coor(1)**2-coor(1) + coor(2)**2-coor(2))
    end function f

    real(rp) function u_ex(x,y)
        real(rp) :: x,y
        u_ex = x*y*(x-1._rp)*(y-1._rp)
    end function u_ex

    !--------------------------------------------!
    !           Transformee affine
    !--------------------------------------------!

    subroutine transformee_triangle(coor_simplexe, coor_triangle)
        ! Subroutine qui permet de transformer les coordonnees d'un triangle simplicial en coordonnees d'un certain triangle
        real(rp), intent(inout), dimension(2,3) :: coor_simplexe ! coordonnees du simplexe qui sont transformees en coordonnees du triangle
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille

        ! Premiere coordonnee
        coor_simplexe(1,:) = coor_triangle(1,1) + (coor_triangle(1,2) - coor_triangle(1,1))*coor_simplexe(1,:)
        coor_simplexe(1,:) = coor_simplexe(1,:) + (coor_triangle(1,3) - coor_triangle(1,1))*coor_simplexe(2,:)
        ! Deuxieme coordonnee
        coor_simplexe(2,:) = coor_triangle(2,1) + (coor_triangle(2,2) - coor_triangle(2,1))*coor_simplexe(1,:) 
        coor_simplexe(2,:) = coor_simplexe(2,:) + (coor_triangle(2,3) - coor_triangle(2,1))*coor_simplexe(2,:)
    end subroutine transformee_triangle

    subroutine init_triangle_ref(coor)
        ! subroutine qui permet d'avoir les coordonnees des sommets du triangle de reference
        real(rp), dimension(2,3), intent(inout) :: coor
        coor(:,1) = (/0._rp, 0._rp/) ! premier sommet
        coor(:,2) = (/1._rp, 0._rp/) ! deuxieme sommet
        coor(:,3) = (/0._rp, 1._rp/) ! troisieme sommet
    end subroutine init_triangle_ref


    real(rp) function det_Jac(coor_triangle)
        ! fonction qui calcule le determinant de la matrice jacobienne de la transformee affine
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        det_Jac = (coor_triangle(1,2)-coor_triangle(1,1))*(coor_triangle(2,3)-coor_triangle(2,1))  
        det_Jac = det_Jac - (coor_triangle(1,3)-coor_triangle(1,1))*(coor_triangle(2,2)-coor_triangle(2,1))
    end function det_Jac


    real(rp) function inv_jac_P1(coor_triangle, k1, k2)
        ! fonction qui calcule le terme (k1,k2) de la matrice inverse de la jacobienne
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        integer :: k1, k2 ! numeros de ligne et colonne du terme de J^-1 qu'on veut recuperer
        real(rp), dimension(2,2) :: inv_jac 
        ! on remplit la matrice
        inv_jac(1,1) = coor_triangle(2,3)-coor_triangle(2,1)
        inv_jac(1,2) = coor_triangle(1,1)-coor_triangle(1,3)
        inv_jac(2,1) = coor_triangle(2,1)-coor_triangle(2,2)
        inv_jac(2,2) = coor_triangle(1,2)-coor_triangle(1,1)
        inv_jac_P1 = inv_jac(k1,k2)*det_Jac(coor_triangle)
    end function inv_jac_P1

    !--------------------------------------------!
    !       Fonctions de base simplexe de ref
    !             et leurs derivees
    !--------------------------------------------!

    real(rp) function phi(coor,i)
        ! fonction qui calcule les fonctions de base P1
        real(rp), dimension(2) :: coor ! coordonnees du point pour lequel on calcule
        integer :: i ! entier qui determine quelle fonction de base on calcule
        ! 1: 1-x-y (sommet S_1(0,0))
        ! 2: x (sommet S_2(1,0))
        ! 3: y (sommet S_3(0,1))
        if (i == 1) then
            phi = 1._rp - coor(1) - coor(2)
        else if (i == 2) then 
            phi = coor(1)
        else if (i == 3) then
            phi = coor(2)
        end if
    end function phi

    real(rp) function dphi_dx(coor,i)
    real(rp), dimension(2) :: coor
    integer :: i
    if (i == 1) then
        dphi_dx = -1._rp
    else if (i == 2) then
        dphi_dx = 1._rp
    else if (i == 3) then
        dphi_dx = 0._rp
    end if
    end function dphi_dx

    real(rp) function dphi_dy(coor,i)
    real(rp), dimension(2) :: coor
    integer :: i
    if (i == 1) then
        dphi_dy = -1._rp
    else if (i == 2) then
        dphi_dy = 0._rp
    else if (i == 3) then
        dphi_dy = 1._rp
    end if
    end function dphi_dy


    !--------------------------------------------!
    !       Formules de quadrature
    !--------------------------------------------!


    subroutine quadrature_triangle_L(quad, coor_triangle)
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        !integer, intent(in) :: num ! numerotation locale du noeud
        real(rp), dimension(2,3) :: coor, coor_ref
        !integer :: j

        call init_triangle_ref(coor)
        coor_ref(:,:) = coor(:,:)

        !on transforme les sommets du triangle de reference en sommets de notre triangle
        call transformee_triangle(coor, coor_triangle) 

        !quad = f(coor(:,1)) * phi(coor_ref(:,1),1) + f(coor(:,2)) * phi(coor_ref(:,2),2) + f(coor(:,3)) * phi(coor_ref(:,3),3) ! quadrature sans le det de la jacob. et le poids de quadrature
        quad = f(coor(:,1)) + f(coor(:,2)) + f(coor(:,3))
        quad = abs(det_Jac(coor_triangle)) * quad * (1._rp/6._rp) ! on multiplie le tout pour avoir notre integrale approchee
    end subroutine quadrature_triangle_L


    subroutine quadrature_triangle_A(quad, coor_triangle)!, num_i, num_j)
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        !integer, intent(in) :: num_i, num_j ! numerotation locale du noeud
        real(rp), dimension(2,3) :: coor
        real(rp) :: poids
        integer :: i
        real(rp) :: A1,A2,A3,A4
        real(rp) :: j11, j12, j21, j22

        poids = (1._rp/6._rp)*abs(det_Jac(coor_triangle))
        ! on recupere les termes de l'inverse de la matrice jacobienne
        j11 = inv_jac_P1(coor_triangle,1,1)
        j12 = inv_jac_P1(coor_triangle,1,2)
        j21 = inv_jac_P1(coor_triangle,2,1)
        j22 = inv_jac_P1(coor_triangle,2,2)

        ! on recupere notre triangle de reference
        call init_triangle_ref(coor)

        quad = 0._rp
        do i = 1,3
            A1 = (j11*dphi_dx(coor(:,i),i)+j12*dphi_dx(coor(:,i),i))
            !A2 = (j11*dphi_dx(coor(:,i),num_j)+j12*dphi_dx(coor(:,i),num_j))
            A3 = (j21*dphi_dy(coor(:,i),i)+j22*dphi_dy(coor(:,i),i))
            !A4 = (j21*dphi_dy(coor(:,i),num_j)+j22*dphi_dy(coor(:,i),num_j))
            !quad = quad + A1*A2 + A3*A4
            quad = quad+A1**2+A3**2
        end do
        
        ! on multiplie le tout par le poids
        quad = quad*poids
    end subroutine quadrature_triangle_A

    !--------------------------------------------!
    !           Remplissage de A et L
    !--------------------------------------------!

    subroutine remplissage_L(L, coordonnees, positions, connect, nb_element, nb_triangle, dim_mat, points_int)
        integer, intent(in) :: nb_element, nb_triangle, dim_mat
        real(rp), dimension(dim_mat), intent(inout) :: L
        real(rp), dimension(nb_element,2), intent(in) :: coordonnees
        integer, dimension(nb_element), intent(in) :: positions
        integer, dimension(3,nb_triangle), intent(in) :: connect
        integer, dimension(dim_mat), intent(in) :: points_int
        integer, dimension(3) :: S ! positions des trois sommets
        real(rp), dimension(nb_element) :: L_avec_frontiere
        integer :: i, j, k ! pour boucle do
        real(rp) :: quad
        real(rp), dimension(2,3) :: coor_triangle

        ! On remplit L avec les points sur la frontiere d'abord
        L_avec_frontiere(:) = 0._rp
        do i = 1,nb_triangle
            S(1) = connect(1,i)
            S(2) = connect(2,i)
            S(3) = connect(3,i)
            do j = 1,3
                if (positions(S(j)) /= 1) then
                    do k = 1,2
                        coor_triangle(k,j) = coordonnees(S(j),k)
                    end do
                    call quadrature_triangle_L(quad, coor_triangle)
                    L_avec_frontiere(S(j)) = L_avec_frontiere(S(j)) + quad
                end if
                !write(6,*) L_avec_frontiere(S(j))
            end do
        end do

        ! Puis on enleve les lignes avec les points du bord pour L
        do i = 1,dim_mat
            L(i) = L_avec_frontiere(points_int(i))
       end do
    end subroutine remplissage_L


    subroutine remplissage_A(A, coordonnees, positions, connect, nb_element, nb_triangle, dim_mat, points_int)
        integer, intent(in) :: nb_element, nb_triangle, dim_mat
        real(rp), dimension(dim_mat, dim_mat), intent(inout) :: A
        real(rp), dimension(nb_element,2), intent(in) :: coordonnees
        integer, dimension(nb_element), intent(in) :: positions
        integer, dimension(3,nb_triangle), intent(in) :: connect
        integer, dimension(dim_mat), intent(in) :: points_int
        real(rp), dimension(nb_element, nb_element) :: A_avec_frontiere
        integer :: ni, nj, m ! pour boucle do
        integer :: i, j
        real(rp) :: quad
        real(rp), dimension(2,3) :: coor_triangle

        ! on fait le calcul
        A_avec_frontiere(:,:) = 0._rp
        do m = 1,nb_triangle ! boucle sur le nombre de triangles
            ! on recupere les coordonnees des sommets du triangle dans lequel on travaille
            do i = 1,3
                do j = 1,2
                    coor_triangle(j,i) = coordonnees(connect(i,m),j)
                end do
            end do
            do ni = 1,3
                i = connect(ni,m)
                if (positions(i) /= 1) then
                    do nj = 1,3
                        j = connect(nj,m)
                        if (positions(j) /= 1) then
                            call quadrature_triangle_A(quad, coor_triangle)!, ni, nj)
                            !write(6,*) quad
                            A_avec_frontiere(i,j) = A_avec_frontiere(i,j) + quad
                        end if
                    end do
                end if
            end do
        end do

        ! Puis on enleve les lignes et colonnes avec les points du bord
        do i = 1,dim_mat
            do j = 1,dim_mat
                A(i,j) = A_avec_frontiere(points_int(i),points_int(j))
            end do
        end do
    end subroutine remplissage_A

    !--------------------------------------------!
    !           Post-traitement
    !--------------------------------------------!

    subroutine recup_vect_U(U, L, points_int, nb_element, dim_mat)
        ! subroutine pour recuperer le vecteur U des solutions avec les points interieurs et sur la frontiere
        integer, intent(in) :: nb_element, dim_mat 
        real(rp), dimension(nb_element), intent(inout) :: U
        real(rp), dimension(dim_mat), intent(in) :: L
        integer, dimension(dim_mat), intent(in) :: points_int
        integer :: i

        U(:) = 0._rp
        do i =1,dim_mat
            U(points_int(i)) = L(i)
        end do
    end subroutine recup_vect_U


end module quadrature_P1