module numerics
    ! Module pour la definition des constantes

    use iso_fortran_env
    IMPLICIT NONE

    !--------------------------------------------!
    !               DANS CE MODULE:
    ! - Definitions des constantes
    ! - Definition du type mat_creuse
    !--------------------------------------------!

    integer, parameter :: rp = real64
    real(rp), parameter :: pi = acos(-1._rp)

    type mat_creuse
        integer :: Ncoefmat ! dimension de Tmat et Jvcell
        integer :: nb_element ! nb_element+1 = dimension de Jposi
        real(rp), dimension(:), allocatable :: Tmat
        integer, dimension(:), allocatable :: JvCell
        integer, dimension(:), allocatable :: Jposi
    end type mat_creuse
end module numerics