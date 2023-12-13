module quadrature_P2
    use numerics
    use quadrature_P1
    IMPLICIT NONE

    contains

    !--------------------------------------------!
    !     Fonctions de base simplexe de ref
    !           et leurs derivees
    !--------------------------------------------!

    real(rp) function phi_P2(coor,i)
        ! fonction qui calcule les fonctions de base P2
        real(rp), dimension(2) :: coor ! coordonnees du point pour lequel on calcule
        integer :: i ! entier qui determine quelle fonction de base on calcule
        ! i = 1, 2, 3, 4, 5 ou 6
        if (i == 1) then
            phi_P2 = 2._rp*phi(coor,1)*(phi(coor,1)-0.5_rp)
        else if (i == 2) then 
            phi_P2 = 2._rp*phi(coor,2)*(phi(coor,2)-0.5_rp)
        else if (i == 3) then
            phi_P2 = 2._rp*phi(coor,3)*(phi(coor,3)-0.5_rp)
        else if (i == 4) then
            phi_P2 = 4._rp*phi(coor, 2)*phi(coor, 3)
        else if (i == 5) then
            phi_P2 = 4._rp*phi(coor, 1)*phi(coor, 3)
        else if (i == 6) then
            phi_P2 = 4._rp*phi(coor, 2)*phi(coor, 1)
        end if
    end function phi_P2

    real(rp) function dphi_dx_P2(coor,i)
    real(rp), dimension(2) :: coor
    integer :: i
    if (i == 1) then
        dphi_dx_P2 = -3._rp + 4._rp*coor(2)+4._rp*coor(1)
    else if (i == 2) then
        dphi_dx_P2 = 1._rp-4._rp*coor(1)
    else if (i == 3) then
        dphi_dx_P2 = 0._rp
    else if (i == 4) then
        dphi_dx_P2 = 4._rp*coor(2)
    else if (i == 5) then
        dphi_dx_P2 = -4._rp*coor(2)
    else if (i == 6) then
        dphi_dx_P2 = 4._rp*coor(1) - 4._rp*coor(2) - 8._rp*coor(1)
    end if
    end function dphi_dx_P2

    real(rp) function dphi_dy_P2(coor,i)
    real(rp), dimension(2) :: coor
    integer :: i
    if (i == 1) then
        dphi_dy_P2 = -3._rp + 4._rp*coor(1)+4._rp*coor(2)
    else if (i == 2) then
        dphi_dy_P2 = 0._rp
    else if (i == 3) then
        dphi_dy_P2 = 1._rp-4._rp*coor(2)
    else if (i == 4) then
        dphi_dy_P2 = 4._rp*coor(1)
    else if (i == 5) then
        dphi_dy_P2 = 4._rp*coor(2) - 4._rp*coor(1) - 8._rp*coor(2)
    else if (i == 6) then
        dphi_dy_P2 = - 4._rp*coor(1)
    end if
    end function dphi_dy_P2

    subroutine init_milieux_triref(coor_milieu)
        ! subroutine qui permet d'avoir les coordonnees des milieux des segments du triangle de reference
        real(rp), dimension(2,3), intent(inout) :: coor_milieu
        coor_milieu(:,1) = (/0.5_rp, 0.5_rp/) ! premier milieu
        coor_milieu(:,2) = (/0._rp, 0.5_rp/) ! deuxieme milieu
        coor_milieu(:,3) = (/0.5_rp, 0._rp/) ! troisieme milieu
    end subroutine init_milieux_triref

    !--------------------------------------------!
    !       Formules de quadrature
    !--------------------------------------------!

    subroutine quadrature_triangle_L_P2(quad, coor_triangle)!, num)
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        !integer, intent(in) :: num ! numerotation locale du noeud
        real(rp), dimension(2,3) :: coor, coor_ref
        real(rp), dimension(2,3) :: coor_milieu_ref, coor_milieu
        !integer :: j

        call init_triangle_ref(coor)
        coor_ref(:,:) = coor(:,:)

        call init_milieux_triref(coor_milieu)
        coor_milieu_ref(:,:) = coor_milieu(:,:)

        !on transforme les sommets du triangle de reference en sommets de notre triangle
        call transformee_triangle(coor, coor_triangle) 
        call transformee_triangle(coor_milieu, coor_triangle)

        quad = f(coor(:,1)) * phi(coor_ref(:,1),1) + f(coor(:,2)) * phi(coor_ref(:,2),2) + f(coor(:,3)) * phi(coor_ref(:,3),3) ! sommets
        quad = quad + f(coor_milieu(:,1)) * phi(coor_milieu_ref(:,1),4) 
        quad = quad + f(coor_milieu(:,2)) * phi(coor_milieu_ref(:,2),5) 
        quad = quad + f(coor_milieu(:,3)) * phi(coor_milieu_ref(:,3),6) 
        quad = abs(det_Jac(coor_triangle)) * quad * (1._rp/12._rp) ! on multiplie le tout pour avoir notre integrale approchee
    end subroutine quadrature_triangle_L_P2


    subroutine quadrature_triangle_A_P2(quad, coor_triangle, num_i, num_j)
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        integer, intent(in) :: num_i, num_j ! numerotation locale du noeud
        real(rp), dimension(2,3) :: coor
        real(rp), dimension(2,3) :: coor_milieu
        real(rp) :: poids
        integer :: i
        real(rp) :: A1,A2,A3,A4
        real(rp) :: j11, j12, j21, j22

        poids = (1._rp/12._rp)*abs(det_Jac(coor_triangle))
        ! on recupere les termes de l'inverse de la matrice jacobienne
        j11 = inv_jac(coor_triangle,1,1)
        j12 = inv_jac(coor_triangle,1,2)
        j21 = inv_jac(coor_triangle,2,1)
        j22 = inv_jac(coor_triangle,2,2)

        ! on recupere notre triangle de reference
        call init_triangle_ref(coor)
        ! on recupere aussi les points milieux
        call init_milieux_triref(coor_milieu)

        quad = 0._rp
        do i = 1,3 ! pour les sommets
            A1 = (j11*dphi_dx_P2(coor(:,i),num_i)+j12*dphi_dx_P2(coor(:,i),num_i))
            A2 = (j11*dphi_dx_P2(coor(:,i),num_j)+j12*dphi_dx_P2(coor(:,i),num_j))
            A3 = (j21*dphi_dy_P2(coor(:,i),num_i)+j22*dphi_dy_P2(coor(:,i),num_i))
            A4 = (j21*dphi_dy_P2(coor(:,i),num_j)+j22*dphi_dy_P2(coor(:,i),num_j))
            quad = quad + A1*A2 + A3*A4
            A1 = (j11*dphi_dx_P2(coor_milieu(:,i),num_i+3)+j12*dphi_dx_P2(coor_milieu(:,i),num_i+3))
            A3 = (j21*dphi_dy_P2(coor_milieu(:,i),num_j+3)+j22*dphi_dy_P2(coor_milieu(:,i),num_j+3))
            quad = quad + A1**2 + A3**2
        end do ! puis pour les points milieux
        !do i = 1,3
        !    A1 = (j11*dphi_dx_P2(coor_milieu(:,i),num_i)+j12*dphi_dx_P2(coor_milieu(:,i),num_i))
        !    A2 = (j11*dphi_dx_P2(coor_milieu(:,i),num_j)+j12*dphi_dx_P2(coor_milieu(:,i),num_j))
        !    A3 = (j21*dphi_dy_P2(coor_milieu(:,i),num_i)+j22*dphi_dy_P2(coor_milieu(:,i),num_i))
        !    A4 = (j21*dphi_dy_P2(coor_milieu(:,i),num_j)+j22*dphi_dy_P2(coor_milieu(:,i),num_j))
        !    quad = quad + A1*A2 + A3*A4
        !end do
        
        ! on multiplie le tout par le poids
        quad = quad*poids
    end subroutine quadrature_triangle_A_P2

end module quadrature_P2