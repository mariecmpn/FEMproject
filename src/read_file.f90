module read_file
    use numerics
    IMPLICIT NONE

    contains

    subroutine open_file(file_name, meshfile_name, save_file, stock, degre, cond)
        character(len=*), intent(in) :: file_name
        character(len=70), intent(inout) :: meshfile_name
        character(len=70), intent(inout) :: save_file
        integer, intent(inout) :: stock
        integer, intent(inout) :: degre
        real(rp), intent(inout) :: cond
        integer :: my_unit = 60

        open(my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')

        read(my_unit,*) meshfile_name
        read(my_unit,*) save_file
        read(my_unit,*) stock
        read(my_unit,*) degre
        read(my_unit,*) cond
    end subroutine open_file

end module read_file