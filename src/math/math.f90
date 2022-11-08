!> d
module math
  use num_types
  use array_math
  use device_math
  use comm
  implicit none

  interface rzero
     module procedure array_rzero_r1, array_rzero_r2, array_rzero_r3, array_rzero_r4
  end interface rzero

contains

  subroutine array_rzero_r1(a, n, bcknd)
    integer, intent(in) :: n    
    real(kind=rp), intent(inout) :: a(:)
    integer, optional :: bcknd
    call array_rzero(a, n)
  end subroutine array_rzero_r1

  subroutine array_rzero_r2(a, n, bcknd)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(:,:)
    integer, optional :: bcknd
    call array_rzero(a, n)
  end subroutine array_rzero_r2

  subroutine array_rzero_r3(a, n, bcknd)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(:,:,:)
    integer, optional :: bcknd
    call array_rzero(a, n)
  end subroutine array_rzero_r3

  subroutine array_rzero_r4(a, n, bcknd)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(:,:,:,:)
    integer, optional :: bcknd
    call array_rzero(a, n)
  end subroutine array_rzero_r4

end module math
