module to_vtk
   use iso_fortran_env
   use numerics
   IMPLICIT NONE
   !integer, parameter :: rp = real64

   type Donnees
      !real(rp) :: Impre
      real(rp), dimension(:), allocatable :: Z
   end type Donnees


   type Meshdef
      integer :: NPoint
      integer :: Nelemt
      !real(rp), dimension(2,1:Npoint) :: coor
      !integer,  dimension(3,1:Nelemt) :: Nu
      real(rp), dimension(:,:), allocatable :: coor
      integer,  dimension(:,:), allocatable :: Nu
   end type Meshdef


   type Variables
      !real(rp), dimension(1:3,:), allocatable :: Ua
      real(rp), dimension(:,:), allocatable :: Ua
   end type Variables   

   contains


   subroutine recup_points(nomfic, nb_pt, nb_tri)
      CHARACTER(len=50) :: nomfic, nomficnode, str, nomficele
      INTEGER :: nb_pt, lens, lens2, nb_tri
      
      lens   = INDEX(nomfic,' ')-1   
      str=''
      lens2  = INDEX(str,' ') 
  
      nomficnode = str
      nomficele = str
      nomficnode(lens2:lens2+lens) = nomfic
      nomficnode(lens2+lens:lens2+lens+7)= '.node'
      nomficele(lens2:lens2+lens) = nomfic
      nomficele(lens2+lens:lens2+lens+6) = '.ele'
  
      OPEN(unit=50,file=nomficnode)
      READ(50,*)nb_pt
      CLOSE(50)
  
      OPEN(unit=50,file=nomficele)
      READ(50,*)nb_tri
      CLOSE(50)
  
  end subroutine recup_points

subroutine recup_nodeel(nomfic, coord, p, t, n_pt, n_tri)

  IMPLICIT NONE

  INTEGER          :: n_node, n_ele, Dim, n1, n2, n3, n_per_triangle
  integer, intent(inout) :: n_pt, n_tri
  INTEGER          :: int, i, lens, lens2!, lens3
  !REAL(rp)          :: r
  INTEGER, DIMENSION(3,n_tri), intent(inout) :: t ! triangles
  INTEGER, DIMENSION(n_pt), intent(inout) :: p ! positions (sur le bord ou non)
  REAL(rp), DIMENSION(n_pt,2), intent(inout) :: coord ! coordonnees
  CHARACTER(len=50), intent(in)        :: nomfic
  CHARACTER(len=50)                   :: nomficnode, nomficele, str


  !--------------------------------------------------------!                                                                                                                                                    
  !    LECTURE DES FICHIERS CARRE.1.NODE ET CARRE.1.ELE    !                                                                                                                                                    
  !--------------------------------------------------------!                                                                                                                                                    
   !print*, 'Nom du fichier .poly'
   !read*,nomfic
  !nomfic = 'carre.3'
  lens   = INDEX(nomfic,' ')-1
  str=''
  lens2  = INDEX(str,' ')

  nomficnode = str ! nom du fichier .node
  nomficele = str ! nom du fichier .ele
  nomficnode(lens2:lens2+lens) = nomfic
  nomficnode(lens2+lens:lens2+lens+7)= '.node'
  nomficele(lens2:lens2+lens) = nomfic
  nomficele(lens2+lens:lens2+lens+6) = '.ele'

  ! lecture des coordonnees et de leur position puis enregistrement dans des tableaux
  OPEN(unit=50,file=nomficnode)
  READ(50,*)n_node,Dim,n1,n2
  DO i=1,n_node
     READ(50,*)int,coord(i,1),coord(i,2),p(i)
  ENDDO

  CLOSE(50)

  OPEN(unit=50,file=nomficele)
  READ(50,*)n_ele,n_per_triangle,n3
  DO i=1,n_ele
     READ(50,*)int,t(1,i),t(2,i),t(3,i)
  ENDDO

  CLOSE(50)
end subroutine recup_nodeel

subroutine points_frontiere (nb_frontiere, n_pt, p, dim_mat, points_int)
   ! On recupere le nombre de points sur la frontiere
   integer, intent(inout) :: nb_frontiere
   integer, intent(in) :: n_pt
   INTEGER, DIMENSION(n_pt), intent(inout) :: p ! positions (sur le bord ou non)
   integer, intent(inout) :: dim_mat
   integer :: i,j
   integer, dimension(:), allocatable, intent(inout) :: points_int

   ! on compte le nombre de points sur la frontiere
   nb_frontiere = 0
   do i = 1,n_pt
      if (p(i) == 1) then
         nb_frontiere = nb_frontiere + 1
      end if
   end do

   ! puis on garde en memoire les numeros des points interieurs
   dim_mat = n_pt - nb_frontiere
   allocate(points_int(dim_mat))
   j = 1
   do i = 1,n_pt
      if (p(i) /= 1) then
         points_int(j) = i
         j = j + 1
      end if
   end do
end subroutine points_frontiere


SUBROUTINE CellVertexVtk_n(DATA, Mesh, PbName)

   TYPE(Donnees)   , INTENT(in)     :: DATA ! 
   TYPE(MeshDef)   , INTENT(in)     :: Mesh !
   !TYPE(Variables) , INTENT(in)     :: Var !
   !CHARACTER(LEN=5),INTENT(in)     :: PbName ! nom du fichier vtk
   CHARACTER(LEN=70), intent(in) :: PbName
   !integer          ,intent(in)     :: flag_topo 

   !REAL(PR)          :: u_x, u_y
   CHARACTER(LEN=75) :: vtk=" "
   INTEGER           :: is, jt
   INTEGER           :: lensd3, lPbName


   WRITE (6,FMT = *) " --------------------------------------------"
   WRITE (6,FMT = *) " ---------      passage dans vtk   ----------"
   WRITE (6,FMT = *) " --------------------------------------------"

   ! ECRITURE sous FICHIER vtk !

   !print*, 'Nom du fichier .vtk'
   !read*,PbName

   lPbName                    = INDEX(PbName,' ') - 1
   vtk(1:lPbName)             = PbName(1:lPbName)
   vtk(lPbName+1:lPbName+5)   = '.vtk '
   lensd3                     = lPbName+4
   OPEN(UNIT=61,FILE=vtk(1:lensd3) )
   
   ! En-tete du fichier vtk
   WRITE(61,'(A)')'# vtk DataFile Version 3.0'
   WRITE(61,'(A)')'# Solution du M1 TR'
   WRITE(61,'(A)')'ASCII'
   WRITE(61,'(A)')'DATASET UNSTRUCTURED_GRID'
   WRITE(61,'(A,I7,A)')'POINTS', Mesh%Npoint,'  float'

   DO is = 1,Mesh%Npoint
      !if (flag_topo == 1) then
      WRITE(61,'(E13.7,2x,E13.7,2x,E13.7)') Mesh%coor(1,is), Mesh%coor(2,is),  DATA%Z(is)
      !else
         !WRITE(61,'(E13.7,2x,E13.7,2x,E13.7)') Mesh%coor(1,is), Mesh%coor(2,is),  Var%Ua(1,is) + DATA%Z(is)
      !end if
   END DO
   WRITE(61,'(A,1x,I7,1x,I6)')'CELLS',Mesh%Nelemt, 4*Mesh%Nelemt

   DO jt = 1,Mesh%Nelemt
      WRITE(61,'(I1,1x,I7,1x,I7,1x,I7)') 3, Mesh%Nu(1,jt)-1, Mesh%Nu(2,jt)-1, Mesh%Nu(3,jt)-1
   END DO

   WRITE(61,'(A,1x,I7)')'CELL_TYPES', Mesh%Nelemt
   DO jt = 1,Mesh%Nelemt
      WRITE(61,'(I1)') 5
   ENDDO

   WRITE(61,'(A,1x,I7)')'POINT_DATA',Mesh%Npoint
   WRITE(61,'(A)')'SCALARS hauteur float'
   WRITE(61,'(A)')'LOOKUP_TABLE DEFAULT'

   !DO is = 1,Mesh%Npoint
   !   WRITE(61,'(ES20.7)') max(1.E-08,Var%Ua(1,is))
   !END DO

   WRITE(61,'(A)')'SCALARS topographie float'
   WRITE(61,'(A)')'LOOKUP_TABLE DEFAULT'

   !DO is = 1,Mesh%Npoint
   !   WRITE(61,'(ES20.7)') max(1.E-08,DATA%Z(is))
   !END DO

   DO is = 1,Mesh%Npoint
      WRITE(61,'(ES20.7)') max(1.E-08,DATA%Z(is))
   END DO

   !WRITE(61,'(A)')'VECTORS vitesse float'

   !DO is = 1,Mesh%Npoint
   !   u_x = vitesse(Var%Ua(1,is) , Var%Ua(2,is))
   !   u_y = vitesse(Var%Ua(1,is) , Var%Ua(3,is))
   !   WRITE(61,'(E13.7,2x,E13.7,2x,E13.7)') max(1.E-08,u_x), max(1.E-08,u_y), 0.
   !END DO

   CLOSE(61)

 END SUBROUTINE CellVertexVtk_n


end module to_vtk