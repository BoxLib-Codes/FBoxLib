module advance_module

  use multifab_module
  use bds_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(phio,phin,dx,dt)

    type(multifab) , intent(in   ) :: phio
    type(multifab) , intent(inout) :: phin
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer i,dm

    ! an array of multifabs; one for each direction
    type(multifab) :: umac(phio%dim)

    ! vector of dx needed for bds routine
    real(kind=dp_t) :: dx_vec(phio%dim)

    logical :: is_conservative(1)

    dm = phio%dim

    dx_vec(:) = dx

    is_conservative(:) = .true.

    do i=1,dm
       call multifab_build_edge(umac(i),phio%la,1,1,i)
       call multifab_setval(umac(i),1.d0,all=.true.)
    end do

    call bds(phio,phin,umac,dx_vec,dt,is_conservative)

    call multifab_fill_boundary(phin)

    ! make sure to destroy the multifab or you'll leak memory
    do i=1,dm
       call multifab_destroy(umac(i))
    end do

  end subroutine advance

end module advance_module

