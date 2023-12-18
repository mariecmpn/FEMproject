module quadrature_P2
    use numerics
    use quadrature_P1
    IMPLICIT NONE

    !--------------------------------------------!
    !               DANS CE MODULE:
    ! - Definition des fonctions de base P2 et 
    ! de leurs derivees
    ! - Formules de quadrature (a 6 points) pour A
    ! et L
    ! - Recherche des points milieux, leur 
    ! position
    ! - Definition de tri_s_m
    !--------------------------------------------!

    contains

    !--------------------------------------------!
    !     Fonctions de base simplexe de ref
    !           et leurs derivees
    !--------------------------------------------!

    real(rp) function phi_P2(coor,i)
        ! fonction qui calcule les fonctions de base P2
        ! S_1(0,0)
        ! S_2(1,0)
        ! S_3(0,1)
        ! S_4(0.5, 0.5)
        ! S_5(0,0.5)
        ! S_6(0.5,0)
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
        ! fonction qui calcule la derivee en x de phi_i, i = 1,..,6
        ! S_1(0,0)
        ! S_2(1,0)
        ! S_3(0,1)
        ! S_4(0.5, 0.5)
        ! S_5(0,0.5)
        ! S_6(0.5,0)
        real(rp), dimension(2) :: coor ! coordonnees du point pour lequel on calcule
        integer :: i ! entier qui determine quelle fonction de base on calcule
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
        ! fonction qui calcule la derivee en y de phi_i, i = 1,..,6
        ! S_1(0,0)
        ! S_2(1,0)
        ! S_3(0,1)
        ! S_4(0.5, 0.5)
        ! S_5(0,0.5)
        ! S_6(0.5,0)
        real(rp), dimension(2) :: coor ! coordonnees du point pour lequel on calcule
        integer :: i ! entier qui determine quelle fonction de base on calcule
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

    subroutine quadrature_triangle_L_P2(quad, coor_triangle, num)
        ! subroutine qui calcule la formule de quadrature a 6 points pour remplir le vecteur L
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        integer, intent(in) :: num ! numerotation locale du noeud
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

        if ((num == 1) .OR. (num == 2) .OR. (num == 3)) then
            quad = f(coor_triangle(:,num))
        else if ((num == 4) .OR. (num == 5) .OR. (num == 6)) then
            quad = f(coor_milieu(:,num-3))
        end if

        !quad = f(coor(:,1)) * phi(coor_ref(:,1),num) + f(coor(:,2)) * phi(coor_ref(:,2),num) + f(coor(:,3)) * phi(coor_ref(:,3),num) ! sommets
        ! puis aux points milieux
        !quad = quad + f(coor_milieu(:,1)) * phi(coor_milieu_ref(:,1),num) 
        !quad = quad + f(coor_milieu(:,2)) * phi(coor_milieu_ref(:,2),num) 
        !quad = quad + f(coor_milieu(:,3)) * phi(coor_milieu_ref(:,3),num) 
        quad = abs(det_Jac(coor_triangle)) * quad * (1._rp/12._rp) ! on multiplie le tout pour avoir notre integrale approchee
    end subroutine quadrature_triangle_L_P2


    subroutine quadrature_triangle_A_P2(quad, coor_triangle, num_i, num_j)
        ! subroutine qui calcule la formule de quadrature a 6 points pour remplir la matrice A
        real(rp), intent(inout) :: quad ! reel qui contiendra la valeur de la quadrature en sortie
        real(rp), intent(in), dimension(2,3) :: coor_triangle ! coordonnees des sommets du triangle dans lequel on travaille
        integer, intent(in) :: num_i, num_j ! numerotation locale du noeud
        real(rp), dimension(2,3) :: coor
        real(rp), dimension(2,3) :: coor_milieu
        real(rp) :: poids
        integer :: i
        real(rp) :: A1,A2,A3,A4
        real(rp) :: j11, j12, j21, j22

        poids = (1._rp/6._rp)*abs(det_Jac(coor_triangle))
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
            A1 = (j11*dphi_dx_P2(coor(:,i),num_i)+j21*dphi_dx_P2(coor(:,i),num_i))
            A2 = (j11*dphi_dx_P2(coor(:,i),num_j)+j21*dphi_dx_P2(coor(:,i),num_j))
            A3 = (j12*dphi_dy_P2(coor(:,i),num_i)+j22*dphi_dy_P2(coor(:,i),num_i))
            A4 = (j12*dphi_dy_P2(coor(:,i),num_j)+j22*dphi_dy_P2(coor(:,i),num_j))
            quad = quad + A1*A2 + A3*A4
        end do ! puis pour les points milieux
        do i = 1,3
            A1 = (j11*dphi_dx_P2(coor_milieu(:,i),num_i)+j21*dphi_dx_P2(coor_milieu(:,i),num_i))
            A2 = (j11*dphi_dx_P2(coor_milieu(:,i),num_j)+j21*dphi_dx_P2(coor_milieu(:,i),num_j))
            A3 = (j12*dphi_dy_P2(coor_milieu(:,i),num_i)+j22*dphi_dy_P2(coor_milieu(:,i),num_i))
            A4 = (j12*dphi_dy_P2(coor_milieu(:,i),num_j)+j22*dphi_dy_P2(coor_milieu(:,i),num_j))
            quad = quad + A1*A2 + A3*A4
        end do
        
        ! on multiplie le tout par le poids
        quad = quad*poids
    end subroutine quadrature_triangle_A_P2


    !--------------------------------------------!
    !       Recherche de points milieux
    !--------------------------------------------!

    subroutine recherche_milieux(coor_milieux, NUBO, Nseg, coordonnees, nb_element)
        ! Subroutine qui permet de trouver les coordonnees milieux de chaque segment 
        integer, intent(in) :: Nseg,nb_element
        integer, dimension(2, Nseg), intent(in) :: NUBO
        real(rp), dimension(Nseg,2), intent(inout) :: coor_milieux
        real(rp), dimension(nb_element,2), intent(in) :: coordonnees
        integer :: i, S1, S2

        coor_milieux(:,:) = 0._rp ! on initialise a 0
        do i = 1,Nseg
            S1 = NUBO(1,i)
            S2 = NUBO(2,i)
            coor_milieux(i,1) = (coordonnees(S1,1) + coordonnees(S2,1))*0.5_rp
            coor_milieux(i,2) = (coordonnees(S1,2) + coordonnees(S2,2))*0.5_rp
        end do 
    end subroutine recherche_milieux

    subroutine pts_et_mid_int(tous_pts_int, dim_matP2, posi_sommet_mid, nb_element, Nseg)
        ! Subroutine qui met dans tous_pts_int les indices des points interieurs (milieux et sommets)
        integer, intent(in) :: nb_element, Nseg
        integer, intent(inout) :: dim_matP2
        integer, dimension(:), allocatable, intent(inout) :: tous_pts_int
        integer, dimension(nb_element+Nseg), intent(in) :: posi_sommet_mid
        integer, dimension(nb_element+Nseg) :: points ! on surestime le nombre de points a l'interieur du domaine
        integer :: i

        dim_matP2 = 0 ! on compte le nombre de points a l'interieur 
        do i = 1,nb_element+Nseg
            if (posi_sommet_mid(i) == 0) then
                dim_matP2 = dim_matP2+1
                points(dim_matP2) = i ! on enregistre dans points l'indice du point interieur 
            end if
        end do
        allocate(tous_pts_int(dim_matP2))
        tous_pts_int(:) = points(1:dim_matP2)
    end subroutine pts_et_mid_int

    subroutine ttes_positions(posi_sommet_mid, positions, nb_element, Nseg, NUBO)
        ! Subroutine qui remplit posi_sommet_mid avec les positions des points (sommets et milieux): 1 si sur le bord, 0 sinon
        integer, dimension(nb_element+Nseg), intent(inout) :: posi_sommet_mid
        integer, dimension(nb_element), intent(in) :: positions
        integer, dimension(2,Nseg) :: NUBO
        integer, intent(in) :: Nseg, nb_element
        integer :: i

        posi_sommet_mid(:) = 0
        posi_sommet_mid(1:nb_element) = positions(:) ! on a deja les positions des sommets

        do i = 1,Nseg ! pour les milieux, on parcourt les segments dans NUBO
            if ((positions(NUBO(1,i)) == 1) .AND. (positions(NUBO(2,i)) == 1)) then 
                posi_sommet_mid(nb_element+i) = 1 ! si les deux sommets du segment sont sur le bord, alors le milieu est aussi sur le bord
            end if 
        end do
    end subroutine ttes_positions


    subroutine connect_milieux(tri_s_m, connect, nb_triangle, NUBO, Nseg, nb_element)
        ! Subroutine qui permet de recuperer les numerotations globales de tous les points (sommets et milieux) d'un triangle
        ! ils sont enregistres dans tri_s_m => de 1 a 3: les numerotations des sommets; de 4 a 6: les numerotations des milieux
        integer, intent(in) :: nb_triangle, Nseg, nb_element
        integer, dimension(6,nb_triangle), intent(inout) :: tri_s_m
        integer, dimension(3,nb_triangle), intent(in) :: connect
        integer, dimension(2,NSeg), intent(in) :: NUBO
        integer :: m, S1, S2, k

        do m = 1,nb_triangle ! boucle sur les triangles
            ! on cherche les indices globaux des sommets et des milieux
            tri_s_m(1:3,m) = connect(:,m) ! on connait deja les sommets grace a connect
            ! premier segment
            if (tri_s_m(1,m) < tri_s_m(2,m)) then
                S1 = tri_s_m(1,m)
                S2 = tri_s_m(2,m)
            else 
                S1 = tri_s_m(2,m)
                S2 = tri_s_m(1,m)
            end if 
            do k = 1,Nseg ! on cherche le segment dans NUBO
                if ((S1 == NUBO(1,k)) .AND. (S2 == NUBO(2,k))) then ! des qu'on le trouve, on note l'indice global du milieu et on sort
                    tri_s_m(6,m) = k+nb_element
                    exit
                end if
            end do
            ! deuxieme segment
            if (tri_s_m(2,m) < tri_s_m(2,m)) then
                S1 = tri_s_m(2,m)
                S2 = tri_s_m(3,m)
            else 
                S1 = tri_s_m(3,m)
                S2 = tri_s_m(2,m)
            end if 
            do k = 1,Nseg ! on cherche le segment dans NUBO
                if ((S1 == NUBO(1,k)) .AND. (S2 == NUBO(2,k))) then ! des qu'on le trouve, on note l'indice global du milieu et on sort
                    tri_s_m(5,m) = k+nb_element
                    exit
                end if
            end do
            ! troisieme segment
            if (tri_s_m(3,m) < tri_s_m(1,m)) then
                S1 = tri_s_m(3,m)
                S2 = tri_s_m(1,m)
            else 
                S1 = tri_s_m(1,m)
                S2 = tri_s_m(3,m)
            end if 
            do k = 1,Nseg ! on cherche le segment dans NUBO
                if ((S1 == NUBO(1,k)) .AND. (S2 == NUBO(2,k))) then ! des qu'on le trouve, on note l'indice global du milieu et on sort
                    tri_s_m(4,m) = k+nb_element
                    exit
                end if
            end do
        end do
    end subroutine connect_milieux

end module quadrature_P2