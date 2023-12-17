module quadrature_P1
    use numerics
    IMPLICIT NONE

    !--------------------------------------------!
    !               DANS CE MODULE:
    ! - Fonctions exactes
    ! - Definition des fonctions de base P1 et 
    ! de leurs derivees
    ! - Definition de la transformee affine, 
    ! le det et l'inverse de sa jacobienne
    ! - Formules de quadrature a 3 points pour
    ! A et L
    !--------------------------------------------!

    contains

    !--------------------------------------------!
    !             Fonctions exactes
    !--------------------------------------------!

    real(rp) function f(coor)
        real(rp), dimension(2) :: coor
        f = 0.5_rp*pi*pi*sin(pi*coor(1))*sin(pi*coor(2))
        !f = 2._rp
    end function f

    real(rp) function u_ex(x,y)
        real(rp) :: x,y
        u_ex = 0.25_rp*sin(pi*x)*sin(pi*y)
    end function u_ex

    !--------------------------------------------!
    !           Transformee affine
    !--------------------------------------------!

    subroutine transformee_triangle(coor_simplexe, coor_triangle)
        ! Subroutine qui permet de transformer les coordonnees d'un triangle simplicial en coordonnees d'un certain triangle
        real(rp), intent(inout), dimension(2,3) :: coor_simplexe ! coordonnees du simplexe qui sont transformees en coordonnees du triangle
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        real(rp), dimension(2,3) :: coor

        coor(:,:) = coor_simplexe(:,:)
        ! Premiere coordonnee
        coor_simplexe(1,:) = coor_triangle(1,1) + (coor_triangle(1,2) - coor_triangle(1,1))*coor(1,:)
        coor_simplexe(1,:) = coor_simplexe(1,:) + (coor_triangle(1,3) - coor_triangle(1,1))*coor(2,:)
        ! Deuxieme coordonnee
        coor_simplexe(2,:) = coor_triangle(2,1) + (coor_triangle(2,2) - coor_triangle(2,1))*coor(1,:) 
        coor_simplexe(2,:) = coor_simplexe(2,:) + (coor_triangle(2,3) - coor_triangle(2,1))*coor(2,:)
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


    real(rp) function inv_jac(coor_triangle, k1, k2)
        ! fonction qui calcule le terme (k1,k2) de la matrice inverse de la jacobienne
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        integer :: k1, k2 ! numeros de ligne et colonne du terme de J^-1 qu'on veut recuperer
        real(rp), dimension(2,2) :: inv_J 
        ! on remplit la matrice
        inv_J(1,1) = coor_triangle(2,3)-coor_triangle(2,1)
        inv_J(1,2) = coor_triangle(1,1)-coor_triangle(1,3)
        inv_J(2,1) = coor_triangle(2,1)-coor_triangle(2,2)
        inv_J(2,2) = coor_triangle(1,2)-coor_triangle(1,1)
        ! pus on renvoit seulement la valeur qui nous interesse 
        inv_jac = inv_J(k1,k2)/det_Jac(coor_triangle)
    end function inv_jac

    real(rp) function mesure_tri(coor_triangle)
        ! fonction qui permet d'obtenir la mesure d'un triangle par la formule de Heron
        real(rp), dimension(2,3) :: coor_triangle
        real(rp) :: demiperi, a, b, c

        a = sqrt((coor_triangle(1,1)- coor_triangle(1,2))**2 + (coor_triangle(2,1)- coor_triangle(2,2))**2)
        b = sqrt((coor_triangle(1,2)- coor_triangle(1,3))**2 + (coor_triangle(2,2)- coor_triangle(2,3))**2)
        c = sqrt((coor_triangle(1,1)- coor_triangle(1,3))**2 + (coor_triangle(2,1)- coor_triangle(2,3))**2)
        demiperi = (a+b+c)*0.5_rp
        mesure_tri = sqrt(demiperi*(demiperi-a)*(demiperi-b)*(demiperi-c))
    end function mesure_tri

    !--------------------------------------------!
    !           Fonctions de reference
    !             et leurs derivees
    !--------------------------------------------!

    real(rp) function phi(coor,i)
        ! fonction qui calcule les fonctions de base P1
        ! 1: 1-x-y (sommet S_1(0,0))
        ! 2: x (sommet S_2(1,0))
        ! 3: y (sommet S_3(0,1))
        real(rp), dimension(2) :: coor ! coordonnees du point pour lequel on calcule
        integer :: i ! entier qui determine quelle fonction de base on calcule
        if (i == 1) then
            phi = 1._rp - coor(1) - coor(2)
        else if (i == 2) then 
            phi = coor(1)
        else if (i == 3) then
            phi = coor(2)
        end if
    end function phi

    real(rp) function dphi_dx(coor,i)
        ! fonction qui calcule la derivee partielle en x de phi_i, i = 1, 2 ou 3
        ! 1: 1-x-y (sommet S_1(0,0))
        ! 2: x (sommet S_2(1,0))
        ! 3: y (sommet S_3(0,1))
        real(rp), dimension(2) :: coor ! coordonnees du point pour lequel on calcule
        integer :: i ! entier qui determine quelle fonction de base on calcule
        if (i == 1) then
            dphi_dx = -1._rp
        else if (i == 2) then
            dphi_dx = 1._rp
        else if (i == 3) then
            dphi_dx = 0._rp
        end if
    end function dphi_dx

    real(rp) function dphi_dy(coor,i)
        ! fonction qui calcule la derivee partielle en y de phi_i, i = 1, 2 ou 3
        ! 1: 1-x-y (sommet S_1(0,0))
        ! 2: x (sommet S_2(1,0))
        ! 3: y (sommet S_3(0,1))
        real(rp), dimension(2) :: coor ! coordonnees du point pour lequel on calcule
        integer :: i ! entier qui determine quelle fonction de base on calcule
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


    subroutine quadrature_triangle_L(quad, coor_triangle, num)
        ! subroutine qui calcule la formule de quadrature a 3 points pour remplir le vecteur L
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        integer, intent(in) :: num ! numerotation locale du noeud
        real(rp), dimension(2,3) :: coor, coor_ref
        !integer :: j

        quad = 0._rp
        !call init_triangle_ref(coor)
        !coor_ref(:,:) = coor(:,:)

        !on transforme les sommets du triangle de reference en sommets de notre triangle
        !call transformee_triangle(coor, coor_triangle) 

        !quad = f(coor(:,1)) * phi(coor_ref(:,1),num) + f(coor(:,2)) * phi(coor_ref(:,2),num) + f(coor(:,3)) * phi(coor_ref(:,3),num) ! quadrature sans le det de la jacob. et le poids de quadrature
        quad = f(coor_triangle(:,num))

        !quad = abs(det_Jac(coor_triangle)) * quad * (1._rp/6._rp)
        quad = quad * mesure_tri(coor_triangle)*(1._rp/3._rp) ! on multiplie le tout pour avoir notre integrale approchee
    end subroutine quadrature_triangle_L


    subroutine quadrature_triangle_A(quad, coor_triangle, num_i, num_j)
        ! subroutine qui calcule la formule de quadrature a 3 points pour remplir la matrice A
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        integer, intent(in) :: num_i, num_j ! numerotation locale du noeud
        real(rp), dimension(2,3) :: coor
        real(rp) :: poids
        real(rp) :: quad1, quad2
        integer :: i, k1
        real(rp) :: A1,A2,A3,A4
        real(rp) :: j11, j12, j21, j22

        ! initialisation des variables
        quad = 0._rp
        quad1 = 0._rp
        quad2 = 0._rp
        poids = (1._rp/2._rp)
        ! on recupere les termes de l'inverse de la matrice jacobienne
        j11 = inv_jac(coor_triangle,1,1)
        j12 = inv_jac(coor_triangle,1,2)
        j21 = inv_jac(coor_triangle,2,1)
        j22 = inv_jac(coor_triangle,2,2)

        ! on recupere notre triangle de reference
        call init_triangle_ref(coor)

        do i = 1,3
            A1 = (j11*dphi_dx(coor(:,i),num_i)+j21*dphi_dx(coor(:,i),num_i))
            A2 = (j11*dphi_dx(coor(:,i),num_j)+j21*dphi_dx(coor(:,i),num_j))
            A3 = (j12*dphi_dy(coor(:,i),num_i)+j22*dphi_dy(coor(:,i),num_i))
            A4 = (j12*dphi_dy(coor(:,i),num_j)+j22*dphi_dy(coor(:,i),num_j))
            quad = quad + A1*A2 + A3*A4
            !write(6,*) 'i = ', i, ' num_i = ', num_i, ' num_j = ', num_j
            !write(6,*) 'dphi_dx(coor(:,i),num_i) = ', dphi_dx(coor(:,i),num_i)
            !write(6,*) 'dphi_dx(coor(:,i),num_j) = ', dphi_dx(coor(:,i),num_j)
            !write(6,*) 'dphi_dy(coor(:,i),num_i) = ', dphi_dy(coor(:,i),num_i)
            !write(6,*) 'dphi_dy(coor(:,i),num_j) = ', dphi_dy(coor(:,i),num_j)
        end do
        !do k1 = 1,2
        !    quad1 = quad1 + inv_jac(coor_triangle,1,k1)*dphi_dx(coor,num_i)
        !    quad2 = quad2 + inv_jac(coor_triangle,1,k1)*dphi_dy(coor,num_j)
        !    quad1 = quad1 + inv_jac(coor_triangle,2,k1)*dphi_dx(coor,num_i)
        !    quad2 = quad2 + inv_jac(coor_triangle,2,k1)*dphi_dy(coor,num_j)
        !    quad = quad + quad1*quad2
        !end do
        
        ! on multiplie le tout par le poids
        quad = quad*poids*abs(det_Jac(coor_triangle))
    end subroutine quadrature_triangle_A




end module quadrature_P1