program test
    use numerics
    use quadrature_P1
    use stockage_matrice
    use to_vtk
    IMPLICIT NONE
    INTEGER :: nb_element, nb_triangle, nb_frontiere, dim_mat
    real(rp), Dimension(:,:), ALLOCATABLE :: coordonnees ! coordonnees des points
    INTEGER, Dimension(:), ALLOCATABLE :: positions ! posiitons des points (1 si sur le bord; 0 sinon)
    integer, DIMENSION(:,:), ALLOCATABLE :: connect ! table de connectivite = table des triangles
    !real(rp), dimension(2,3) :: coor_triangle
    !real(rp) :: quad = 0._rp
    CHARACTER(len=70) :: nom = 'carre.3'
    real(rp), dimension(:), allocatable :: L, U ! vecteurs L second membre et U 
    real(rp), dimension(:, :), allocatable :: A ! matrice de masse A
    integer :: i,j
    integer, dimension(:), allocatable :: points_int ! vecteur avec les numeros des points interieurs
    integer :: Nseg ! nb de segments
    integer, dimension(:,:), allocatable :: NUBO
    real(rp), dimension(:), allocatable :: Tmat
    integer, dimension(:), allocatable :: Jposi
    integer, dimension(:), allocatable :: JvCell
    integer :: Ncoefmat

    !external :: DGETRF, DGETRS ! pour le 
    !external :: DPPTRF, DPPTRS
    !external :: DPOTRF, DPOTRS
    !external :: dgesv

    !external :: DPOSV 
    !integer :: info
    !integer, dimension(:), allocatable :: ipiv
  
    ! recuperation nombre de points pour allouer les tableaux
    call recup_points(nom, nb_element, nb_triangle)

    ! allocation tableaux
    allocate(positions(nb_element), coordonnees(nb_element,2), connect(3,nb_triangle))

    ! recuperation des tableaux dans les fichiers .node et .ele
    call recup_nodeel(nom, coordonnees, positions, connect, nb_element, nb_triangle)

    ! recuperation du nombre de points sur la frontiere
    call points_frontiere(nb_frontiere, nb_element, positions, dim_mat, points_int)

    !do i = 1,nb_triangle
    !    write(6,*) (connect(j,i), j = 1,3)
    !end do

    ! Stockage creux matrice A


    ! allocation de L, A et U
    !allocate(L(dim_mat), A(dim_mat, dim_mat), U(nb_element))

    !------------------------------------------------------------

    call remplissage_NUBO_cher(NUBO, nb_triangle, connect, Nseg)

    write(6,*) 'NUBO : '
    do i = 1,Nseg
        write(6,*) (NUBO(j,i), j = 1,2)
    end do
    

    call calcul_NcoefMat(NUBO, Nseg, NcoefMat, nb_element)

    write(6,*) 'NcoefMat = ', NcoefMat

    allocate(Tmat(1:NcoefMat)) 
    allocate(Jposi(1:nb_element+1)) 
    allocate(JvCell(1:NcoefMat))

    call tableaux_assemblage_mat(Jposi, JvCell, Ncoefmat, nb_element, NUBO, Nseg)

    call stockage_morse(Tmat, Jposi, JvCell, coordonnees, connect, nb_element, nb_triangle, NcoefMat)

    write(6,*) 'Tmat : '
    write(6,*) (Tmat(i), i = 1,NcoefMat)

    write(6,*) 'Jposi : '
    write(6,*) (Jposi(i), i = 1,nb_element+1)

    write(6,*) 'JvCell : '
    write(6,*) (JvCell(i), i = 1,NcoefMat)

    !------------------------------------------------------------

    !call remplissage_NUBO_peu_cher(NUBO, nb_element, nb_triangle, connect, Nseg)

    
    ! remplissage vecteur L
    !call remplissage_L(L, coordonnees, positions, connect, nb_element, nb_triangle, dim_mat, points_int)
    
    !do i = 1,dim_mat
    !    write(6,*) L(i)
    !end do

    ! remplissage matrice A
    !call remplissage_A(A, coordonnees, positions, connect, nb_element, nb_triangle, dim_mat, points_int)

    !do i = 1,dim_mat
        !write(6,*) (A(i,j), j = 1,dim_mat)
    !end do

    !allocate(ipiv(dim_mat))
    !do i = 1,dim_mat
        !ipiv(i) = i
    !end do

    !call DPOTRS('U',dim_mat,1,A,dim_mat,L,dim_mat,info)
    !call DGETRF(dim_mat,dim_mat,A,dim_mat,ipiv,info)
    !call DGETRS('N',dim_mat,1,A,dim_mat,ipiv,L,dim_mat,info)
    !call DPPTRF('U', dim_mat, A, info)
    !call DPPTRS('U', dim_mat, 1, A, L, dim_mat, info)
    !call DPOTRF('U',dim_mat,A,dim_mat,info)
    !call DPOTRS('U',dim_mat,1,A,dim_mat,L,dim_mat,info)
    !call dgesv(dim_mat,1,A,dim_mat,ipiv,L,dim_mat,info)

    !call DPOSV('U',dim_mat,1,A,dim_mat,L,dim_mat,info)

    !write(6,*)
    !do i = 1,dim_mat
        !write(6,*) (A(i,j), j = 1,dim_mat)
    !end do

    !write(6,*) info
    !do i = 1,dim_mat
    !    write(6,*) L(i)
    !end do

    !call recup_vect_U(U, L, points_int, nb_element, dim_mat)

    !write(6,*)
    !do i = 1,nb_element
        !write(6,*) U(i)
    !end do

    deallocate(coordonnees, positions, connect, NUBO, Tmat, Jposi, JvCell)

end program test