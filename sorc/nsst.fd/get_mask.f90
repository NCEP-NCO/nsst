 subroutine get_mask(filename,mask,im,jm)
!
! abstract: read 1/12 degree land/water mask and convert to conventional values
!
! Input/output
 character(*), intent(in) :: filename
 integer     , intent(in) :: im,jm
 integer, dimension(im,jm) :: mask
! Local
 integer, dimension(im,jm) :: mask0
 integer, parameter :: lun_rtgmsk=31


 open(lun_rtgmsk,file=filename,form='formatted')
 read (lun_rtgmsk,31) ((mask0(i,j),i=1,im),j=1,jm)
 31 format (80I1)
!
! reset mask = 0 for ocean; mask = 1 for land
!
  do j=1,jm
     do i=1,im
        if( mask0(i,j) == 0) mask(i,j) = 0
        if( mask0(i,j) == 3) mask(i,j) = 1
     enddo
  enddo

 close (31)

 end subroutine get_mask
