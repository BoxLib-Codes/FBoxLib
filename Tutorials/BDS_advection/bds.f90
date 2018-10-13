module bds_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: bds

contains

  subroutine bds(so,sn,umac,dx,dt,is_conserv)

    type(multifab) , intent(in   ) :: so
    type(multifab) , intent(inout) :: sn
    type(multifab) , intent(in   ) :: umac(:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_conserv(:)

    ! this will hold slx, sly, and slxy
    type(multifab) :: slope

    real(kind=dp_t), pointer ::    sop(:,:,:,:)
    real(kind=dp_t), pointer ::    snp(:,:,:,:)
    real(kind=dp_t), pointer :: slopep(:,:,:,:)
    real(kind=dp_t), pointer ::  uadvp(:,:,:,:)
    real(kind=dp_t), pointer ::  vadvp(:,:,:,:)
    real(kind=dp_t), pointer ::  wadvp(:,:,:,:)

    integer :: dm,ng_s,ng_c,ng_u,n,i,comp
    integer :: lo(so%dim),hi(so%dim)

    dm = so%dim

    if (dm .eq. 2) then
       ! 3 components and 1 ghost cell
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slxy
       call multifab_build(slope,so%la,3,1)
    else if (dm .eq. 3) then
       ! 7 components and 1 ghost cell
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slz
       ! component 4 = slxy
       ! component 5 = slxz
       ! component 6 = slyz
       ! component 7 = slxyz
       call multifab_build(slope,so%la,7,1)
    end if

    ng_s = so%ng
    ng_c = slope%ng
    ng_u = umac(1)%ng

    do i=1,nfabs(so)
       sop    => dataptr(so, i)
       snp    => dataptr(sn, i)
       slopep => dataptr(slope, i)
       uadvp  => dataptr(umac(1), i)
       vadvp  => dataptr(umac(2), i)
       lo =  lwb(get_box(so, i))
       hi =  upb(get_box(so, i))
       select case (dm)
       case (2)
          ! only advancing the tracer
          do comp=1,so%nc
             call bdsslope_2d(lo, hi, &
                              sop(:,:,1,comp), ng_s, &
                              slopep(:,:,1,:), ng_c, &
                              dx) 

             call bdsconc_2d(lo, hi, &
                             sop(:,:,1,comp), snp(:,:,1,comp), ng_s, &
                             slopep(:,:,1,:), ng_c, &
                             uadvp(:,:,1,1), vadvp(:,:,1,1), ng_u, &
                             dx, dt, is_conserv(comp))
          end do
       case (3)
          wadvp  => dataptr(umac(3), i)
          ! only advancing the tracer
          do comp=1,so%nc
             call bdsslope_3d(lo, hi, &
                              sop(:,:,:,comp), ng_s, &
                              slopep(:,:,:,:), ng_c, &
                              dx)

             call bdsconc_3d(lo, hi, &
                             sop(:,:,:,comp), snp(:,:,:,comp), ng_s, &
                             slopep(:,:,:,:), ng_c, &
                             uadvp(:,:,:,1), vadvp(:,:,:,1), wadvp(:,:,:,1), ng_u, &
                             dx, dt, is_conserv(comp))
          end do
       end select
    end do

    call multifab_destroy(slope)

  end subroutine bds

  subroutine bdsslope_2d(lo,hi,s,ng_s,slope,ng_c,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local variables
    real(kind=dp_t), allocatable :: sint(:,:)

    real(kind=dp_t) :: diff(4)
    real(kind=dp_t) :: smin(4)
    real(kind=dp_t) :: smax(4)
    real(kind=dp_t) :: sc(4)

    real(kind=dp_t) :: hx,hy,eps
    real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
    integer         :: i,j,ll,mm

    logical :: limit_slopes

    limit_slopes = .false.

    ! nodal with one ghost cell
    allocate(sint(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2))

    hx = dx(1)
    hy = dx(2)

    eps = 1.d-10

    ! bicubic interpolation to corner points
    ! (i,j,k) refers to lower corner of cell
    do j = lo(2)-1,hi(2)+2
       do i = lo(1)-1,hi(1)+2
          sint(i,j) = (s(i-2,j-2) + s(i-2,j+1) + s(i+1,j-2) + s(i+1,j+1) &
               - 7.d0*(s(i-2,j-1) + s(i-2,j  ) + s(i-1,j-2) + s(i  ,j-2) + & 
                       s(i-1,j+1) + s(i  ,j+1) + s(i+1,j-1) + s(i+1,j  )) &
              + 49.d0*(s(i-1,j-1) + s(i  ,j-1) + s(i-1,j  ) + s(i  ,j  )) ) / 144.d0
       enddo
    enddo

    do j = lo(2)-1,hi(2)+1
       do i = lo(1)-1,hi(1)+1 

          ! compute initial estimates of slopes from unlimited corner points

          ! sx
          slope(i,j,1) = 0.5d0*(sint(i+1,j+1) + sint(i+1,j  ) - &
                                sint(i  ,j+1) - sint(i  ,j  ) ) / hx

          ! sy
          slope(i,j,2) = 0.5d0*(sint(i+1,j+1) - sint(i+1,j  ) + &
                                sint(i  ,j+1) - sint(i  ,j  ) ) / hy

          ! sxy
          slope(i,j,3) = ( sint(i+1,j+1) - sint(i+1,j  ) &
                          -sint(i  ,j+1) + sint(i  ,j  ) ) / (hx*hy)
          
          if (limit_slopes) then

             ! ++ / sint(i+1,j+1)
             sc(4) = s(i,j) + 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2))  &
                  + 0.25d0*hx*hy*slope(i,j,3)

             ! +- / sint(i+1,j  )
             sc(3) = s(i,j) + 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2))  &
                  - 0.25d0*hx*hy*slope(i,j,3)

             ! -+ / sint(i  ,j+1)
             sc(2) = s(i,j) - 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2)) &
                  - 0.25d0*hx*hy*slope(i,j,3)

             ! -- / sint(i  ,j  )
             sc(1) = s(i,j) - 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2)) &
                  + 0.25d0*hx*hy*slope(i,j,3)

             ! enforce max/min bounds
             smin(4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smax(4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             
             smin(3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smax(3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             
             smin(2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smax(2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             
             smin(1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
             smax(1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
             
             do mm=1,4
                sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
             enddo

             ! iterative loop
             do ll = 1,3 
                sumloc = 0.25d0*(sc(4) + sc(3) + sc(2) + sc(1))
                sumdif = (sumloc - s(i,j))*4.d0
                sgndif = sign(1.d0,sumdif)

                do mm=1,4
                   diff(mm) = (sc(mm) - s(i,j))*sgndif
                enddo

                kdp = 0

                do mm=1,4
                   if (diff(mm) .gt. eps) then
                      kdp = kdp+1
                   end if
                end do

                do mm = 1,4 
                   if (kdp.lt.1) then 
                      div = 1.d0
                   else
                      div = dble(kdp)
                   end if

                   if (diff(mm).gt.eps) then
                      redfac = sumdif*sgndif/div
                      kdp = kdp-1
                   else
                      redfac = 0.d0
                   end if

                   if (sgndif .gt. 0.d0) then
                      redmax = sc(mm) - smin(mm)
                   else
                      redmax = smax(mm) - sc(mm)
                   end if

                   redfac = min(redfac,redmax)
                   sumdif = sumdif - redfac*sgndif
                   sc(mm) = sc(mm) - redfac*sgndif
                enddo
             enddo

             ! final slopes

             ! sx
             slope(i,j,1) = 0.5d0*( sc(4) + sc(3) &
                  -sc(1) - sc(2))/hx

             ! sy
             slope(i,j,2) = 0.5d0*( sc(4) + sc(2) &
                  -sc(1) - sc(3))/hy

             ! sxy
             slope(i,j,3) = ( sc(1) + sc(4) &
                  -sc(2) - sc(3) ) / (hx*hy)

          end if

       enddo
    enddo

    deallocate(sint)

  end subroutine bdsslope_2d

  subroutine bdsslope_3d(lo,hi,s,ng_s,slope,ng_c,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_c
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(inout) :: slope(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local variables
    real(kind=dp_t), allocatable :: sint(:,:,:)

    real(kind=dp_t) :: diff(8)
    real(kind=dp_t) :: smin(8)
    real(kind=dp_t) :: smax(8)
    real(kind=dp_t) :: sc(8)

    real(kind=dp_t) :: c1,c2,c3,c4
    real(kind=dp_t) :: hx,hy,hz,eps
    real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
    integer         :: i,j,k,ll,mm

    logical :: limit_slopes

    limit_slopes = .false.

    ! nodal with one ghost cell
    allocate(sint(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2))

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    eps = 1.d-10

    c1 = (343.d0/1728.d0)
    c2 = (49.d0 /1728.d0)
    c3 = (7.d0  /1728.d0)
    c4 = (1.d0  /1728.d0)

    ! tricubic interpolation to corner points
    ! (i,j,k) refers to lower corner of cell
    do k = lo(3)-1,hi(3)+2
       do j = lo(2)-1,hi(2)+2
          do i = lo(1)-1,hi(1)+2
             sint(i,j,k) = c1*( s(i  ,j  ,k  ) + s(i-1,j  ,k  ) + s(i  ,j-1,k  ) &
                               +s(i  ,j  ,k-1) + s(i-1,j-1,k  ) + s(i-1,j  ,k-1) &
                               +s(i  ,j-1,k-1) + s(i-1,j-1,k-1) ) &
                          -c2*( s(i-1,j  ,k+1) + s(i  ,j  ,k+1) + s(i-1,j-1,k+1) &
                               +s(i  ,j-1,k+1) + s(i-1,j+1,k  ) + s(i  ,j+1,k  ) &
                               +s(i-2,j  ,k  ) + s(i+1,j  ,k  ) + s(i-2,j-1,k  ) &
                               +s(i+1,j-1,k  ) + s(i-1,j-2,k  ) + s(i  ,j-2,k  ) &
                               +s(i-1,j+1,k-1) + s(i  ,j+1,k-1) + s(i-2,j  ,k-1) &
                               +s(i+1,j  ,k-1) + s(i-2,j-1,k-1) + s(i+1,j-1,k-1) &
                               +s(i-1,j-2,k-1) + s(i  ,j-2,k-1) + s(i-1,j  ,k-2) &
                               +s(i  ,j  ,k-2) + s(i-1,j-1,k-2) + s(i  ,j-1,k-2) ) &
                          +c3*( s(i-1,j+1,k+1) + s(i  ,j+1,k+1) + s(i-2,j  ,k+1) &
                               +s(i+1,j  ,k+1) + s(i-2,j-1,k+1) + s(i+1,j-1,k+1) &
                               +s(i-1,j-2,k+1) + s(i  ,j-2,k+1) + s(i-2,j+1,k  ) &
                               +s(i+1,j+1,k  ) + s(i-2,j-2,k  ) + s(i+1,j-2,k  ) &
                               +s(i-2,j+1,k-1) + s(i+1,j+1,k-1) + s(i-2,j-2,k-1) &
                               +s(i+1,j-2,k-1) + s(i-1,j+1,k-2) + s(i  ,j+1,k-2) &
                               +s(i-2,j  ,k-2) + s(i+1,j  ,k-2) + s(i-2,j-1,k-2) &
                               +s(i+1,j-1,k-2) + s(i-1,j-2,k-2) + s(i  ,j-2,k-2) ) &
                          -c4*( s(i-2,j+1,k+1) + s(i+1,j+1,k+1) + s(i-2,j-2,k+1) &
                               +s(i+1,j-2,k+1) + s(i-2,j+1,k-2) + s(i+1,j+1,k-2) &
                               +s(i-2,j-2,k-2) + s(i+1,j-2,k-2) )
          enddo
       enddo
    enddo

    do k = lo(3)-1,hi(3)+1
       do j = lo(2)-1,hi(2)+1
          do i = lo(1)-1,hi(1)+1 

             ! compute initial estimates of slopes from unlimited corner points

             ! sx
             slope(i,j,k,1) = 0.25d0*( ( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                                        +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                                        +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) / hx

             ! sy
             slope(i,j,k,2) = 0.25d0*( ( sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) &
                                        +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                        +sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)) ) / hy

             ! sz
             slope(i,j,k,3) = 0.25d0*( ( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                                        +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                        +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) / hz

             ! sxy
             slope(i,j,k,4) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j  ,k+1) &
                                       +sint(i+1,j+1,k  ) + sint(i+1,j+1,k+1)) &
                                     -( sint(i+1,j  ,k  ) + sint(i+1,j  ,k+1) &
                                       +sint(i  ,j+1,k  ) + sint(i  ,j+1,k+1)) ) / (hx*hy)

             ! sxz
             slope(i,j,k,5) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                                       +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                                     -( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                                       +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) / (hx*hz)

             ! syz
             slope(i,j,k,6) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                       +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                     -( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                                       +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) / (hy*hz)

             ! sxyz
             slope(i,j,k,7) = (-sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) + sint(i  ,j+1,k  ) &
                               +sint(i  ,j  ,k+1) - sint(i+1,j+1,k  ) - sint(i+1,j  ,k+1) &
                               -sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) ) / (hx*hy*hz)

             if (limit_slopes) then

             ! +++ / sint(i+1,j+1,k+1)
             sc(8) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! ++- / sint(i+1,j+1,k  )
             sc(7) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! +-+ / sint(i+1,j  ,k+1)
             sc(6) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! +-- / sint(i+1,j  ,k  )
             sc(5) = s(i,j,k) &
                  +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! -++ / sint(i  ,j+1,k+1)
             sc(4) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! -+- / sint(i  ,j+1,k  )
             sc(3) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! --+ / sint(i  ,j  ,k+1)
             sc(2) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                  +0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! ---/ sint(i  ,j  ,k  )
             sc(1) = s(i,j,k) &
                  +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                  +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                  -0.125d0*hx*hy*hz*slope(i,j,k,7)

             ! enforce max/min bounds
             smin(8) = min(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                           s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
             smax(8) = max(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                           s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))

             smin(7) = min(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                           s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))
             smax(7) = max(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                           s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))

             smin(6) = min(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                           s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))
             smax(6) = max(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                           s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))

             smin(5) = min(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                           s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))
             smax(5) = max(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                           s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))

             smin(4) = min(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                           s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))
             smax(4) = max(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                           s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))

             smin(3) = min(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                           s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))
             smax(3) = max(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                           s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))

             smin(2) = min(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                           s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))
             smax(2) = max(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                           s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))

             smin(1) = min(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                           s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))
             smax(1) = max(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                           s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))

             do mm=1,8
                sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
             enddo
             
             ! iterative loop
             do ll = 1,3 
                sumloc = 0.125d0*(sc(1)+sc(2)+sc(3)+sc(4)+sc(5)+sc(6)+sc(7)+sc(8))
                sumdif = (sumloc - s(i,j,k))*8.d0
                sgndif = sign(1.d0,sumdif)

                do mm=1,8
                   diff(mm) = (sc(mm) - s(i,j,k))*sgndif
                enddo

                kdp = 0

                do mm=1,8
                   if (diff(mm) .gt. eps) then
                      kdp = kdp+1
                   end if
                end do

                do mm = 1,8
                   if (kdp.lt.1) then 
                      div = 1.d0
                   else
                      div = dble(kdp)
                   end if

                   if (diff(mm).gt.eps) then
                      redfac = sumdif*sgndif/div
                      kdp = kdp-1
                   else
                      redfac = 0.d0
                   end if

                   if (sgndif .gt. 0.d0) then
                      redmax = sc(mm) - smin(mm)
                   else
                      redmax = smax(mm) - sc(mm)
                   end if

                   redfac = min(redfac,redmax)
                   sumdif = sumdif - redfac*sgndif
                   sc(mm) = sc(mm) - redfac*sgndif
                enddo
             enddo

             ! final slopes

             ! sx
             slope(i,j,k,1) = 0.25d0*( ( sc(5) + sc(7) &
                                        +sc(6) + sc(8)) &
                                      -( sc(1) + sc(3) &
                                        +sc(2) + sc(4)) ) / hx

             ! sy
             slope(i,j,k,2) = 0.25d0*( ( sc(3) + sc(7) &
                                        +sc(4) + sc(8)) &
                                      -( sc(1) + sc(5) &
                                        +sc(2) + sc(6)) ) / hy

             ! sz
             slope(i,j,k,3) = 0.25d0*( ( sc(2) + sc(6) &
                                        +sc(4) + sc(8)) &
                                      -( sc(1) + sc(5) &
                                        +sc(3) + sc(7)) ) / hz

             ! sxy
             slope(i,j,k,4) = 0.5d0*( ( sc(1) + sc(2) &
                                       +sc(7) + sc(8)) &
                                     -( sc(5) + sc(6) &
                                       +sc(3) + sc(4)) ) / (hx*hy)

             ! sxz
             slope(i,j,k,5) = 0.5d0*( ( sc(1) + sc(3) &
                                       +sc(6) + sc(8)) &
                                     -( sc(5) + sc(7) &
                                       +sc(2) + sc(4)) ) / (hx*hz)

             ! syz
             slope(i,j,k,6) = 0.5d0*( ( sc(1) + sc(5) &
                                       +sc(4) + sc(8)) &
                                     -( sc(2) + sc(6) &
                                       +sc(3) + sc(7)) ) / (hy*hz)

             ! sxyz
             slope(i,j,k,7) = (-sc(1) + sc(5) + sc(3) &
                               +sc(2) - sc(7) - sc(6) &
                               -sc(4) + sc(8) ) / (hx*hy*hz)

             endif

          enddo
       enddo
    enddo


  end subroutine bdsslope_3d

  subroutine bdsconc_2d(lo,hi,s,sn,ng_s,slope,ng_c,uadv,vadv,ng_u,dx,dt,is_conserv)

    integer        ,intent(in   ) :: lo(:),hi(:),ng_s,ng_c,ng_u
    real(kind=dp_t),intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t),intent(inout) ::    sn(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t),intent(in   ) :: slope(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t),intent(in   ) ::  uadv(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t),intent(in   ) ::  vadv(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t),intent(in   ) :: dx(:),dt
    logical        ,intent(in   ) :: is_conserv

    ! local variables
    real(kind=dp_t),allocatable ::   siphj(:,:)
    real(kind=dp_t),allocatable ::   sijph(:,:)

    real(kind=dp_t) :: hx,hy,hxs,hys,gamp,gamm
    real(kind=dp_t) :: vtrans,stem,vaddif,vdif
    real(kind=dp_t) :: isign, jsign
    real(kind=dp_t) :: uconv,vconv,divu
    real(kind=dp_t) :: u1,u2,v1,v2,uu,vv

    integer i,j,iup,jup

    allocate(siphj(lo(1):hi(1)+1,lo(2):hi(2)  ))
    allocate(sijph(lo(1):hi(1)  ,lo(2):hi(2)+1))

    hx = dx(1)
    hy = dx(2)

    do j = lo(2),hi(2) 
       do i = lo(1)-1,hi(1) 

          ! ******************************* 
          ! calculate Gamma plus for flux F

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j+1)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j
             jsign = 1.d0
             u2 = uadv(i+1,j)
          else 
             jup   = j+1
             jsign = -1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j+1) .gt. 0) then
                u2 = uadv(i+1,j+1)
             end if
          end if

          vv = vadv(iup,j+1)

          hxs = hx*isign
          hys = hy*jsign

          gamp = s(iup,jup)+     &
               (hxs*.5 - (u1+u2)*dt/3.d0)*slope(iup,jup,1) +   &
               (hys*.5 -    vv*dt/3.d0)*slope(iup,jup,2) +    &
               (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
               vv*(2.*u2+u1)*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma plus for flux F
          ! ****************************************

          ! *****************************************
          ! calculate Gamma minus for flux F

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j-1
             jsign = 1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j-1) .gt. 0) then
                u2 = uadv(i+1,j-1)
             end if
          else 
             jup   = j
             jsign = -1.d0
             u2 = uadv(i+1,j)
          end if

          vv = vadv(iup,j)

          hxs = hx*isign
          hys = hy*jsign

          gamm = s(iup,jup)+     &
               (hxs*.5 - (u1+u2)*dt/3.d0)*slope(iup,jup,1) +   &
               (hys*.5 -    vv*dt/3.d0)*slope(iup,jup,2) +    &
               (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
               vv*(2.*u2+u1)*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma minus for flux F
          ! ****************************************

          ! *********************************
          ! calculate siphj

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vdif = 0.5d0*dt*(vadv(iup,j+1)*gamp -  &
               vadv(iup,j)*gamm ) / hy
          stem = s(iup,j) + (isign*hx - uadv(i+1,j)*dt)*0.5d0*slope(iup,j,1)
          vaddif = stem*0.5d0*dt*( &
               uadv(iup+1,j) - uadv(iup,j))/hx
          divu = (uadv(iup+1,j)-uadv(iup,j))/hx +  &
                 (vadv(iup,j+1)-vadv(iup,j))/hy 
          siphj(i+1,j) = stem - vdif - vaddif + 0.5d0*dt*stem*divu

          ! end of calculation of siphj
          ! *************************************

       enddo
    enddo


    do j = lo(2)-1,hi(2) 
       do i = lo(1),hi(1)

          ! ********************************** 
          ! calculate Gamma plus for flux G

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i+1,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i
             isign = 1.d0
             v2 = vadv(i,j+1)
          else
             iup   = i+1
             isign = -1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i+1,j+1) .gt. 0) then
                v2 = vadv(i+1,j+1)
             end if
          end if

          uu = uadv(i+1,jup)       

          hxs = hx*isign
          hys = hy*jsign

          gamp = s(iup,jup)+ &
               (hys*.5 - (v1+v2)*dt/3.)*slope(iup,jup,2) +   &
               (hxs*.5 - uu*dt/3.)*slope(iup,jup,1) + &
               (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
               (2.*v2+v1)*uu*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma plus for flux G
          ! ****************************************

          ! *****************************************
          ! calculate Gamma minus for flux G

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i-1
             isign = 1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i-1,j+1) .gt. 0) then
                v2 = vadv(i-1,j+1)
             end if
          else
             iup   = i
             isign = -1.d0
             v2 = vadv(i,j+1)
          end if

          uu = uadv(i,jup)       

          hxs = hx*isign
          hys = hy*jsign

          gamm = s(iup,jup) +    &
               (hys*.5 - (v1+v2)*dt/3.)*slope(iup,jup,2) +    &
               (hxs*.5 - uu*dt/3.)*slope(iup,jup,1) +   &
               (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
               (2.*v2+v1)*uu*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma minus for flux G
          ! ****************************************

          ! *********************************
          ! calculate sijph

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vdif = 0.5d0*dt* &
               (uadv(i+1,jup)*gamp-uadv(i,jup)*gamm)/hx
          stem = s(i,jup) + (jsign*hy - vadv(i,j+1)*dt)*0.5d0*slope(i,jup,2)
          vaddif = stem*0.5d0*dt*(vadv(i,jup+1) - vadv(i,jup))/hy
          divu =  (uadv(i+1,jup)-uadv(i,jup))/hx +  &
               (vadv(i,jup+1)-vadv(i,jup))/hy 
          sijph(i,j+1) = stem - vdif - vaddif + 0.5d0*dt*stem*divu

          ! end of calculation of sijph
          ! *************************************

       enddo
    enddo

    ! advance solution
    if (is_conserv) then

       ! conservative update
       do j = lo(2),hi(2) 
          do i = lo(1),hi(1) 
             sn(i,j) = s(i,j) - dt*(  &
                  (siphj(i+1,j)*uadv(i+1,j)-siphj(i,j)*uadv(i,j))/hx +  &
                  (sijph(i,j+1)*vadv(i,j+1)-sijph(i,j)*vadv(i,j))/hy)
          enddo
       enddo

    else 

       ! non-conservative update
       do j = lo(2),hi(2) 
          do i = lo(1),hi(1) 
             uconv = 0.5d0 * (uadv(i+1,j)+uadv(i,j))
             vconv = 0.5d0 * (vadv(i,j+1)+vadv(i,j))

             sn(i,j) = s(i,j) - dt*( &
                  uconv * (siphj(i+1,j) - siphj(i,j)) / hx + &
                  vconv * (sijph(i,j+1) - sijph(i,j)) / hy )
          enddo
       enddo

    endif

    deallocate(siphj,sijph)

  end subroutine bdsconc_2d

  subroutine bdsconc_3d(lo,hi,s,sn,ng_s,slope,ng_c,uadv,vadv,wadv,ng_u,dx,dt,is_conserv)

    integer        ,intent(in   ) :: lo(:),hi(:),ng_s,ng_c,ng_u
    real(kind=dp_t),intent(in   ) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t),intent(inout) ::    sn(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t),intent(in   ) :: slope(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t),intent(in   ) ::  uadv(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t),intent(in   ) ::  vadv(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t),intent(in   ) ::  wadv(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t),intent(in   ) :: dx(:),dt
    logical        ,intent(in   ) :: is_conserv

    ! local variables
    integer i,j,k,ioff,joff,koff,ll

    real(kind=dp_t), allocatable :: sedgex(:,:,:)
    real(kind=dp_t), allocatable :: sedgey(:,:,:)
    real(kind=dp_t), allocatable :: sedgez(:,:,:)

    real(kind=dp_t), allocatable :: ux(:,:,:)
    real(kind=dp_t), allocatable :: vy(:,:,:)
    real(kind=dp_t), allocatable :: wz(:,:,:)

    real(kind=dp_t) :: isign,jsign,ksign,hx,hy,hz,uconv,vconv,wconv
    real(kind=dp_t) :: del(3),p1(3),p2(3),p3(3),p4(3)
    real(kind=dp_t) :: val1,val2,val3,val4,val5
    real(kind=dp_t) :: u,v,w,uu,vv,ww,gamma,gamma2
    real(kind=dp_t) :: dt2,dt3,dt4,half,sixth

    allocate(sedgex(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ))
    allocate(sedgey(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ))
    allocate(sedgez(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1))

    allocate(ux(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(wz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    dt2 = dt/2.d0
    dt3 = dt/3.d0
    dt4 = dt/4.d0

    half = 0.5d0
    sixth = 1.d0/6.d0

    ! compute cell-centered ux, vy, and wz
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             ux(i,j,k) = (uadv(i+1,j,k) - uadv(i,j,k)) / hx
             vy(i,j,k) = (vadv(i,j+1,k) - vadv(i,j,k)) / hy
             wz(i,j,k) = (wadv(i,j,k+1) - wadv(i,j,k)) / hz
          end do
       end do
    end do

    ! compute sedgex on x-faces
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgex without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j,k) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ! centroid of rectangular volume
             del(1) = isign*0.5d0*hx - 0.5d0*uadv(i,j,k)*dt
             del(2) = 0.d0
             del(3) = 0.d0
             call eval(s(i+ioff,j,k),slope(i+ioff,j,k,:),del,sedgex(i,j,k))

             ! source term
             sedgex(i,j,k) = sedgex(i,j,k) - dt2*sedgex(i,j,k)*ux(i+ioff,j,k)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k) .gt. 0) then
                u = uadv(i,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j+joff,k) + gamma*vy(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j+1,k)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j+1,k)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vadv(i+ioff,j+1,k)
             sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k) .gt. 0) then
                u = uadv(i,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j+joff,k) + gamma*vy(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vadv(i+ioff,j,k)
             sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j,k+koff) .gt. 0) then
                u = uadv(i,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k+1)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j,k+koff) + gamma*wz(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k+1)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k+1)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i+ioff,j,k+1)
             sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j,k+koff) .gt. 0) then
                u = uadv(i,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j,k+koff) + gamma*wz(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i+ioff,j,k)
             sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hz)

          enddo
       enddo
    enddo

    ! compute sedgey on y-faces    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgey without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             ! centroid of rectangular volume
             if (vadv(i,j,k) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             del(1) = 0.d0
             del(2) = jsign*0.5d0*hy - 0.5d0*vadv(i,j,k)*dt
             del(3) = 0.d0
             call eval(s(i,j+joff,k),slope(i,j+joff,k,:),del,sedgey(i,j,k))

             ! source term
             sedgey(i,j,k) = sedgey(i,j,k) - dt2*sedgey(i,j,k)*vy(i,j+joff,k)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k) .gt. 0) then
                v = vadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i+ioff,j+joff,k) + gamma*ux(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i+1,j+joff,k)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i+1,j+joff,k)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
             gamma = gamma * uadv(i+1,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k) .gt. 0) then
                v = vadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i+ioff,j+joff,k) + gamma*ux(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i,j+joff,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i,j+joff,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * uadv(i,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i,j,k+koff) .gt. 0) then
                v = vadv(i,j,k+koff)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i,j+joff,k+koff) + gamma*wz(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k+1)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k+1)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i,j+joff,k+1)
             sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i,j,k+koff) .gt. 0) then
                v = vadv(i,j,k+koff)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i,j+joff,k+koff) + gamma*wz(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hz)

          enddo
       enddo
    enddo

    ! compute sedgez on z-faces
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgez without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             ! centroid of rectangular volume
             if (wadv(i,j,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             del(1) = 0.d0
             del(2) = 0.d0
             del(3) = ksign*0.5d0*hz - 0.5d0*wadv(i,j,k)*dt
             call eval(s(i,j,k+koff),slope(i,j,k+koff,:),del,sedgez(i,j,k))

             ! source term
             sedgez(i,j,k) = sedgez(i,j,k) - dt2*sedgez(i,j,k)*wz(i,j,k+koff)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j,k) .gt. 0) then
                w = wadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = isign*0.5d0*hx - uadv(i+1,j,k+koff)*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i+ioff,j,k+koff) + gamma*ux(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i+1,j,k+koff)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i+1,j,k+koff)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * uadv(i+1,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j,k) .gt. 0) then
                w = wadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = isign*0.5d0*hx - uadv(i,j,k+koff)*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i+ioff,j,k+koff) + gamma*ux(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i,j,k+koff)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i,j,k+koff)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * uadv(i,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i,j+joff,k) .gt. 0) then
                w = wadv(i,j+joff,k)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - vadv(i,j+1,k+koff)*dt
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i,j+joff,k+koff) + gamma*vy(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j+1,k+koff)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j+1,k+koff)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
             gamma = gamma * vadv(i,j+1,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i,j+joff,k) .gt. 0) then
                w = wadv(i,j+joff,k)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k+koff)*dt
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i,j+joff,k+koff) + gamma*vy(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j,k+koff)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j,k+koff)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vadv(i,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hy)

          enddo
       enddo
    enddo

    ! advance solution
    if (is_conserv) then

       ! conservative update
       do k = lo(3),hi(3)
          do j = lo(2),hi(2) 
             do i = lo(1),hi(1) 
                sn(i,j,k) = s(i,j,k) - dt*(  &
                     (sedgex(i+1,j,k)*uadv(i+1,j,k)-sedgex(i,j,k)*uadv(i,j,k))/hx +  &
                     (sedgey(i,j+1,k)*vadv(i,j+1,k)-sedgey(i,j,k)*vadv(i,j,k))/hy + &
                     (sedgez(i,j,k+1)*wadv(i,j,k+1)-sedgez(i,j,k)*wadv(i,j,k))/hz )
             enddo
          enddo
       enddo

    else 

       ! non-conservative update
       do k = lo(3),hi(3)
          do j = lo(2),hi(2) 
             do i = lo(1),hi(1) 
                uconv = 0.5d0 * (uadv(i+1,j,k)+uadv(i,j,k))
                vconv = 0.5d0 * (vadv(i,j+1,k)+vadv(i,j,k))
                wconv = 0.5d0 * (wadv(i,j,k+1)+wadv(i,j,k))

                sn(i,j,k) = s(i,j,k) - dt*( &
                     uconv * (sedgex(i+1,j,k) - sedgex(i,j,k)) / hx + &
                     vconv * (sedgey(i,j+1,k) - sedgey(i,j,k)) / hy + &
                     wconv * (sedgez(i,j,k+1) - sedgez(i,j,k)) / hz )
             enddo
          enddo
       enddo

    endif

    deallocate(sedgex,sedgey,sedgez)
    deallocate(ux,vy,wz)

  end subroutine bdsconc_3d

  subroutine eval(s,slope,del,val)

    real(kind=dp_t), intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: slope(:)
    real(kind=dp_t), intent(in   ) :: del(:)
    real(kind=dp_t), intent(  out) :: val

    val = s + del(1)*slope(1) + del(2)*slope(2) + del(3)*slope(3) &
         + del(1)*del(2)*slope(4) + del(1)*del(3)*slope(5) + del(2)*del(3)*slope(6) &
         + del(1)*del(2)*del(3)*slope(7)

  end subroutine eval

end module bds_module
