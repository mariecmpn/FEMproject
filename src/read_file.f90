module read_file
    use numerics
    IMPLICIT NONE

    contains

    subroutine open_file(file_name, meshfile_name, save_file, stock, degre, cond, res)
        ! Subroutine qui recupere les informations du fichier file_name
        character(len=*), intent(in) :: file_name ! nom du fichier a lire
        character(len=70), intent(inout) :: meshfile_name ! nom du fichier de maillage qu'on utilisera
        character(len=70), intent(inout) :: save_file ! nom du fichier vtk dans lequel on enregistrera la solution
        integer, intent(inout) :: stock ! type de stockage (en P1): creux ou normal
        integer, intent(inout) :: degre ! degre des elements finis: P1 ou P2
        real(rp), intent(inout) :: cond ! condition de Dirichlet (constante)
        character(len = 1), intent(inout) :: res ! type de resolution: gradient conjugue ou lapack
        integer :: my_unit = 60

        ! ouverture du fichier
        open(my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')

        ! Lecture
        read(my_unit,*) meshfile_name
        read(my_unit,*) save_file
        read(my_unit,*) stock
        read(my_unit,*) degre
        read(my_unit,*) cond
        read(my_unit,*) res

        ! Fermeture
        close(my_unit)
    end subroutine open_file

end module read_file