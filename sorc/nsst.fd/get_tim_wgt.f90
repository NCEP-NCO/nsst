subroutine get_tim_wei(iy,im,id,ih,mon1,mon2,wei1,wei2)
!$$$  
! Abstract:  get two months and weights for monthly climatology
! Created by Xu Li, March, 2019

 implicit none

! Input
 integer, intent(in) :: iy,im,id,ih
! Output
 integer, intent(out) :: mon1,mon2
 real,    intent(out) :: wei1,wei2

! Local declare
 real :: rjday
 integer :: mon,monend,monm,monp,jdow,jdoy,jday
 integer :: jda(8)
!
!dayhf : julian day of the middle of each month
!
 real, dimension(13) ::  dayhf
 data dayhf/15.5,45.0,74.5,105.0,135.5,166.0,196.5,227.5,258.0,288.5,319.0,349.5,380.5/

! 15, 44, 73.5, 104, 134.5, 165, 195.5, 226.5, 257, 287.5, 318.5, 349  ! from WOA05

 jda=0
 jda(1)=iy
 jda(2)=im
 jda(3)=id
 jda(5)=ih
 jdow = 0
 jdoy = 0
 jday = 0
 call w3doxdat(jda,jdow,jdoy,jday)
 rjday=jdoy+jda(5)/24.
 if(rjday.lt.dayhf(1)) rjday=rjday+365.

 monend = 12
 do mon = 1, monend
    monm = mon
    monp = mon + 1
    if( rjday >= dayhf(monm) .and. rjday < dayhf(monp) ) then
       mon1 = monm
       mon2 = monp
       go to 10
    endif
 enddo

 print *,'wrong rjday',rjday
 call abort
 10     continue

 wei1 = (dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
 wei2 = (rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))

 if( mon2 == 13 ) mon2=1

 write(*,'(a,2I4,3F9.3)') 'mon1,mon2,rjday,wei1,wei2=',mon1,mon2,rjday,wei1,wei2

 end subroutine get_tim_wei
