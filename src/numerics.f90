module numerics
    ! Module pour la definition des constantes

    use iso_fortran_env
    IMPLICIT NONE

    integer, parameter :: rp = real64
    real(rp), parameter :: pi = acos(-1._rp)
end module numerics