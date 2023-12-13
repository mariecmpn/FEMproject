module stockage_matrice
    use numerics
    use quadrature_P1
    use quadrature_P2
    use remplissage_A_L
    IMPLICIT NONE

    ! Declaration d’un type maillon
    Type maillon
        integer :: valeur ! Numero du sommet
        type(maillon), Pointer :: suivant ! Pointeur vers le maillon suivant 
    End Type maillon

    contains

    !--------------------------------------------!
    !   Remplissage vecteur NUBO des segments
    !               du maillage
    !--------------------------------------------!

    subroutine remplissage_NUBO_cher(NUBO, nb_triangle, connect, Nseg)
        ! remplissage du tableau NUBO avec l'algorithme simple mais cher
        integer, intent(in) :: nb_triangle
        integer, dimension(:,:), allocatable, intent(inout) :: NUBO ! tableau qui contient les sommets de tous les segments
        integer, dimension(:,:), allocatable :: NUBO_entier
        integer, dimension(3,nb_triangle), intent(in) :: connect
        integer, intent(inout) :: Nseg ! nombre de segments
        integer :: t, i, j, k, s1, s2, s, j1
        logical :: pas_segment

        ! on surestime le nombre de segments
        Nseg = 3*nb_triangle
        allocate(NUBO_entier(2,Nseg))
        NUBO_entier(:,:) = 0

        t = 0
        do i = 1,nb_triangle ! on parcourt les triangles
            do j = 1,3 ! on parcourt les sommets du triangle
                ! on recupere les indices globaux des sommets j et j+1 (mod 3)
                s1 = connect(j,i)
                if (j == 3) then
                    j1 = 1
                else 
                    j1 = j+1
                end if
                s2 = connect(j1,i)
                ! on ordonne les indices
                s = 0
                if (s1 > s2) then
                    s = s1
                    s1 = s2
                    s2 = s
                end if
                if (t == 0) then 
                    ! si premier segment on le met direct
                    NUBO_entier(1,t+1) = s1
                    NUBO_entier(2,t+1) = s2
                    t = t+1
                else ! sinon on parcourt le tableau NUBO pour voir s'il y est deja
                    k = 1
                    pas_segment = .TRUE.
                    do while (pas_segment .AND. (k<=t))
                        if (s1 == NUBO_entier(1,k) .AND. s2 == NUBO_entier(2,k)) then ! il y est deja, on sort de la boucle
                            pas_segment = .FALSE.
                        else ! il n'y est pas, on continue a parcourir le tableau
                            k = k+1
                        end if
                    end do
                    if (pas_segment) then ! il n'est pas deja dans le tableau, on le rajoute
                        NUBO_entier(1,t+1) = s1
                        NUBO_entier(2,t+1) = s2
                        t = t+1
                    end if
                end if
            end do
        end do

        ! enfin on enleve les lignes vides
        Nseg = t
        allocate(NUBO(2,Nseg))
        NUBO(:,:) = 0
        NUBO(:,:) = NUBO_entier(:,1:t)
    end subroutine remplissage_NUBO_cher


    subroutine remplissage_NUBO_peu_cher(NUBO, nb_element, nb_triangle, connect, Nseg)
        ! Remplissage de NUBO avec l'algorithme moins cher mais plus complique
        integer, dimension(:,:), allocatable, intent(inout) :: NUBO ! tableau qui contient les sommets de tous les segments
        integer, intent(in) :: nb_element, nb_triangle
        integer, dimension(3,nb_triangle), intent(in) :: connect
        integer, intent(inout) :: Nseg
        integer, dimension(nb_element) :: NSupVois
        Type(maillon), dimension(:), allocatable :: Hashtab 
        integer :: i, j ! pour boucles do
        integer :: ismin, ismax, s, t, j1
        type(maillon), pointer :: next
        type(maillon) :: new_maillon
        logical :: sortie
        
        allocate(Hashtab(1:nb_element))

        do i=1,nb_element
            NULLIFY(Hashtab(i)%suivant) ! Tous les maillons de Hashtab pointent sur le pointeur NULL
            Hashtab(i)%valeur=i
        end do

        ! Initialisations de Nseg et NSupVois 
        Nseg=0
        NSupVois(1:nb_element)=0

        do i = 1,nb_triangle ! on parcourt les triangles
            do j = 1,3 ! on parcourt les sommets des triangles
                ! on recupere les indices globaux des sommets j et j+1 (mod 3)
                ismin = connect(j,i)
                j1 = j+1
                if (j1 == 4) j1 = modulo(j,3)
                ismax = connect(j1,i)
                ! on ordonne les indices
                if (ismax < ismin) then
                    ismin = s
                    ismin = ismax
                    ismax = s
                end if
                ! on parcourt la liste chainee Hashtab(ismin)
                next = Hashtab(ismin)%suivant
                sortie = .TRUE.
                do while (sortie)
                    if (next%valeur == ismax) then ! le maillon existe deja: on sort de la boucle
                        sortie = .FALSE.
                    else if (next%valeur > ismax) then ! le numero du maillon suivant est > ismax: on l'insere dans la liste
                        new_maillon%valeur = ismax
                        new_maillon%suivant = next
                        next%suivant = new_maillon
                        sortie = .FALSE.
                        ! on rajoute un voisin
                        Nseg=Nseg+1
                        NSupVois(ismin)=NSupVois(ismin)+1
                    else if (.NOT. associated(next)) then ! on arrive en bout de chaine sans etre tombe sur ismax: on rajoute le maillon a la fin
                        new_maillon%valeur = ismax
                        new_maillon%suivant = next
                        next%suivant = new_maillon
                        sortie = .FALSE.
                        ! on rajoute un voisin
                        Nseg=Nseg+1
                        NSupVois(ismin)=NSupVois(ismin)+1
                    else ! sinon on continue a parcourir la liste
                        next = next%suivant 
                    end if
                end do
            end do
        end do 

        allocate(NUBO(2,Nseg))

        t=0
        do i=1,nb_element ! Boucle sur tous les sommets (numerotes globalement)
            if (NSupVois(i) /= 0) then ! Si le sommet i a des voisins
                ! On recupere les voisins superieurs en parcourant la liste chaınee Hashtab(i) 
                next = Hashtab(ismin)%suivant
                do s = 1,NSupVois(i)
                    j = next%valeur
                    ! et on remplit NUBO au fur et a mesure
                    t=t+1
                    NUBO(1,t)=i
                    NUBO(2,t)=j ! j est le numero du maillon sur lequel on est
                    next = next%suivant 
                end do
            end if 
        end do
    end subroutine remplissage_NUBO_peu_cher

    !--------------------------------------------!
    !     Procedures pour le stockage Morse
    !               de la matrice A
    !--------------------------------------------!

    subroutine calcul_NcoefMat(NUBO, Nseg, NcoefMat, nb_element)
        ! Subroutine pour le calcul du coefficient NCoefMat
        integer, intent(in) :: nb_element, Nseg
        integer, intent(inout) :: NcoefMat
        integer, dimension(2, Nseg), intent(in) :: NUBO
        integer :: iseg, js, is

        ! Calcul de NcoefMat
        NcoefMat=nb_element ! contribution des termes diagonaux a_ii 
        do iseg=1,Nseg
            is=NUBO(1,iseg) 
            js=NUBO(2,iseg)
            if (is<=nb_element .and. js<=nb_element) then
                NcoefMat=NcoefMat+2
            end if
        end do
    end subroutine calcul_NCoefMat


    subroutine tableaux_assemblage_mat(Jposi, JvCell, Ncoefmat, nb_element, NUBO, Nseg)
        ! Subroutine qui remplit les tableaux Jposi et JvCell
        integer, intent(in) :: nb_element, Nseg
        integer, intent(in) :: NcoefMat
        integer, dimension(2, Nseg), intent(in) :: NUBO
        !real(rp), dimension(1:NcoefMat), intent(inout) :: Tmat
        integer, dimension(1:NcoefMat), intent(inout) :: JvCell
        integer, dimension(1:nb_element+1), intent(inout) :: Jposi
        integer, dimension(:), allocatable :: IndPL
        integer :: iseg, js, is, jv, kv, tmp

        ! Allocation des tableaux 
        !allocate(Tmat(1:NcoefMat)) 
        !allocate(Jposi(1:nb_element+1)) 
        !allocate(JvCell(1:NcoefMat))

        do is=1,nb_element+1
            Jposi(is)=is
        end do
        do iseg=1,Nseg
            is=NUBO(1,iseg)
            js=NUBO(2,iseg)
            if (is<=nb_element .and. js<=nb_element) then ! si on est sur un segment interieur
                Jposi(is+1:nb_element+1)=Jposi(is+1:nb_element+1)+1
                Jposi(js+1:nb_element+1)=Jposi(js+1:nb_element+1)+1
            end if
        end do

        ! Initialisation de IndPL
        allocate(IndPL(1:nb_element))
        IndPL(1:nb_element)=Jposi(1:nb_element)

        ! Indices de colonnes des termes diagonaux
        do is=1,nb_element
            JvCell(IndPL(is))=is
            IndPL(is)=IndPL(is)+1
        end do

        ! Indices de colonnes des termes extra-diagonaux
        do iseg=1,Nseg
            is=NUBO(1,iseg) 
            js=NUBO(2,iseg)
            if (is<=nb_element .and. js<=nb_element) then
                JvCell(IndPL(is))=js
                IndPL(is)=IndPL(is)+1
                JvCell(IndPL(js))=is
                IndPL(js)=IndPL(js)+1
            end if 
        end do

        do is=1,nb_element ! Boucle sur les lignes
            do jv=Jposi(is+1)-1,Jposi(is),-1 ! Boucle descendante sur les elements non nuls de la ligne
                do kv=Jposi(is)+1,jv ! Boucle sur les premiers termes
                    if (JvCell(kv-1)>JvCell(kv)) then 
                        tmp=JvCell(kv-1) 
                        JvCell(kv-1)=JvCell(kv) 
                        JvCell(kv)=tmp
                    end if 
                end do
            end do 
        end do
    end subroutine tableaux_assemblage_mat


    function AmatLoc(coor_triangle)
        ! Fonction qui calcule la matrice AmatLoc pour un triangle donne
        real(rp), dimension(2,3), intent(in) :: coor_triangle
        real(rp), dimension(3,3) :: AmatLoc
        integer :: ni, nj

        do ni=1,3 ! Boucle sur les 3 sommets numerotes localement
            do nj=1,3 ! idem
                ! Calcul de a_ij^m
                call quadrature_triangle_A(AmatLoc(ni,nj), coor_triangle, ni, nj)
            end do 
        end do
    end function AmatLoc


    subroutine stockage_morse(Tmat, Jposi, JvCell, coordonnees, connect, nb_element, nb_triangle, NcoefMat)
        ! Subroutine pour le remplissage de Tmat
        integer, intent(in) :: nb_element, nb_triangle, NcoefMat
        real(rp), dimension(1:NcoefMat), intent(inout) :: Tmat
        integer, dimension(1:NcoefMat), intent(inout) :: JvCell
        integer, dimension(1:nb_element+1), intent(inout) :: Jposi
        real(rp), dimension(nb_element, 2), intent(in) :: coordonnees
        integer, dimension(3,nb_triangle), intent(in) :: connect
        real(rp), dimension(3,3) :: Aloc
        real(rp), dimension(2,3) :: coor_triangle
        integer :: i,j,m

        
        do m=1,nb_triangle
            do i = 1,3
                do j = 1,2
                    coor_triangle(j,i) = coordonnees(connect(i,m),j)
                end do
            end do
            Aloc = AmatLoc(coor_triangle)
            call Assemble(m, connect, Aloc, Tmat, Jposi, JvCell, nb_element, nb_triangle, NcoefMat)
        end do

    end subroutine stockage_morse


    subroutine Assemble(m, connect, Aloc, Tmat, Jposi, JvCell, nb_element, nb_triangle, NcoefMat)
        ! Subroutine qui ajoute la contribution Aloc a Tmat aux bons endroits
        integer, intent(in) :: m, nb_element, nb_triangle, NcoefMat
        integer, dimension(3,nb_triangle), intent(in) :: connect
        real(rp), dimension(3,3), intent(in) :: Aloc
        real(rp), dimension(1:NcoefMat), intent(inout) :: Tmat
        integer, dimension(1:NcoefMat), intent(in) :: JvCell
        integer, dimension(1:nb_element+1), intent(in) :: Jposi
        integer :: i,j,k
        ! On recupere les numeros globaux des sommets du triangle 
        i = connect(1,m)
        j = connect(2,m)
        k = connect(3,m)
        ! On ajoute les contributions ligne par ligne ! Pour la ligne i, on a
        if (i<=nb_element) then ! on ne touche pas le bord
            ! on ajoute dans A la contribution du triangle au terme a(phi_i,phi_i) :
            call Ajout(i,i,Aloc(1,1),Tmat,Jposi,JvCell, nb_element, NcoefMat) 
        end if
        if (j<=nb_element) then ! on ne touche pas le bord
            ! on ajoute dans A la contribution du triangle au terme a(phi_i,phi_j) :
            call Ajout(i,j,Aloc(1,2),Tmat,Jposi,JvCell,nb_element, NcoefMat) 
        end if
        if (k<=nb_element) then ! on ne touche pas le bord
            ! on ajoute dans A la contribution du triangle au terme a(phi_i,phi_k) : 
            call Ajout(i,k,Aloc(1,3),Tmat,Jposi,JvCell, nb_element, NcoefMat)
        end if
        ! Remarque : si on doit ajouter des termes de bord, on le fait ici avec un else end if

        ! Pour la ligne j, on recommence une procedure similaire
        if (i<=nb_element) then ! on ne touche pas le bord
            call Ajout(j,i,Aloc(2,1),Tmat,Jposi,JvCell, nb_element, NcoefMat) 
        end if
        if (j<=nb_element) then ! on ne touche pas le bord
            call Ajout(j,j,Aloc(2,2),Tmat,Jposi,JvCell, nb_element, NcoefMat) 
        end if
        if (k<=nb_element) then ! on ne touche pas le bord
            call Ajout(j,k,Aloc(2,3),Tmat,Jposi,JvCell, nb_element, NcoefMat)
        end if

        ! Pour la ligne k, on recommence une procedure similaire
        if (i<=nb_element) then ! on ne touche pas le bord
            call Ajout(k,i,Aloc(3,1),Tmat,Jposi,JvCell, nb_element, NcoefMat) 
        end if
        if (j<=nb_element) then ! on ne touche pas le bord
            call Ajout(k,j,Aloc(3,2),Tmat,Jposi,JvCell, nb_element, NcoefMat) 
        end if
        if (k<=nb_element) then ! on ne touche pas le bord
            call Ajout(k,k,Aloc(3,3),Tmat,Jposi,JvCell, nb_element, NcoefMat)
        end if
    end subroutine Assemble


    subroutine Ajout(ii, jj, coefA, Tmat, Jposi, JvCell, nb_element, NcoefMat)
        ! Subroutine qui ajoute, pour les indices ii et jj, coefA dans Tmat 
        integer, intent(in) :: nb_element, NcoefMat
        integer, intent(in) :: ii, jj
        real(rp), intent(in) :: coefA
        real(rp), dimension(1:NcoefMat), intent(inout) :: Tmat
        integer, dimension(1:NcoefMat), intent(in) :: JvCell
        integer, dimension(1:nb_element+1), intent(in) :: Jposi
        integer :: j
        logical :: trouver

        trouver=.false.
        do j = Jposi(ii),Jposi(ii+1)-1 
            if (JvCell(j) == jj) then
                Tmat(j) = Tmat(j)+coefA
                trouver = .true.
                exit ! on a ajoute notre contribution donc on peut arreter de chercher ou la mettre
            end if
        end do
        if (trouver.eqv..false.) then
            print*,"Probleme d’assemblage de la matrice A"
            stop 
        end if
    end subroutine Ajout

    !--------------------------------------------!
    !     Reconstruction de A grace a son
    !               stockage Morse
    !--------------------------------------------!

    subroutine reconstruction_A_morse(A, Tmat, Jposi, JvCell, points_int, nb_element, NcoefMat, dim_mat)
        integer, intent(in) :: nb_element, NcoefMat, dim_mat
        real(rp), dimension(1:NcoefMat), intent(in) :: Tmat
        integer, dimension(1:NcoefMat), intent(in) :: JvCell
        integer, dimension(1:nb_element+1), intent(in) :: Jposi
        real(rp), dimension(dim_mat,dim_mat), intent(inout) :: A
        real(rp), dimension(nb_element,nb_element) :: A_avec_frontiere
        integer, dimension(1:dim_mat), intent(in) :: points_int
        integer :: i,j,jj,nb_non0

        A_avec_frontiere(:,:) = 0._rp
        jj = 1 
        do i = 1,nb_element
            ! On regarde le nombre de termes non nuls pour la ligne i
            nb_non0 = Jposi(i+1) - Jposi(i)
            !print*, nb_non0
            do j = 1,nb_non0 
                A_avec_frontiere(i,JvCell(jj)) = Tmat(jj) ! puis on remplit au bon endroit pour chaque terme non nul
                jj = jj+1
            end do
        end do

        ! Enfin on ne garde que les points interieurs
        do i = 1,dim_mat
            do j = 1,dim_mat
                A(i,j) = A_avec_frontiere(points_int(i),points_int(j))
            end do
        end do
    end subroutine reconstruction_A_morse

end module stockage_matrice