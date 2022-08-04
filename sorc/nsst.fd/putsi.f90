 subroutine putsi(anl,mask,nx,ny)
!
! abstract: Fill land grids with derived temperature
! 
! author: Xu Li, March, 2019, based on putsi in RTG SST analysis
!
! input:
!       anl      : temperature analysis
!       mask     : surface mask
!       nx       : dimension in x-direction
!       ny       : dimension in y-direction
! output:
!       anl      : temperature analysis, include land grids

 real,    dimension(nx,ny), intent(inout) :: anl
 integer, dimension(nx,ny), intent(inout) :: mask
 integer,                   intent(in   ) :: nx,ny
! Local:
 real :: weight,dlon,dlat
 integer :: landflag

 do j = 1, ny
    do i = 1, nx
       if( mask(i,j) /= 1 ) then
          if ( anl(i,j) < -1.8 ) anl(i,j) = -1.8
       endif
    enddo
 enddo

 landflag = 1
 weight = 1.45
 dlon = 360.0 / real(nx)
 dlat = 180.0 / real(ny)
!
 call sor_lapl(anl, mask, nx, ny, dlat, dlon, landflag, weight)

 end subroutine putsi
