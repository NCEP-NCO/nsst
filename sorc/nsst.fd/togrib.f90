 subroutine togrib (fname,anal,mask,maskflag,imx,imy,iyear,imon,iday,ihr)
! abstract: writes out the dates, oi analysis and oi error in grib
!
! author: Xu Li (based on togrib in RTG SST analysis)  ORG: W/NP20     date: 2019-03-13
!
! program history log:
!
! usage:  call togrib(fname,anal,mask,maskflag,imx,imy,iyear,imon,iday)
!
!
!   INPUT ARGUMENT LIST:
!     fname    - file name to store the result in grib
!     anal     - array of sst analysis field (deg K)
!     mask     - array containing land/sea mask
!     maskflag - 1 = with bitmap; 0 = w/o bit map
!     imx      - first dimension of anal array
!     imy      - second dimension of anal array
!     iyear    - year (4 digits) of analysis
!     imon     - month of analysis
!     iday     - day of analysis
!     ihr      - hour of analysis
!
! attributes:
!   language: fortran 90
!
 implicit none
 character(*), intent(in) :: fname
 real,    dimension(imx,imy), intent(in) :: anal
 integer, dimension(imx,imy), intent(in) :: mask
 integer, intent(in) :: maskflag,imx,imy,iyear,imon,iday,ihr
! Local
 integer, dimension(200)      :: kpds=0
 integer, dimension(200)      :: kgds=0
 real,    dimension(imx,imy) :: out,out_mask
 logical*1, dimension(imx,imy) :: lb
 integer :: icen,iyr,kf,i,j,jj,iret
 integer :: latinc,loninc
 integer, parameter :: lun_grb=11

 kf=imx*imy

 iyr    = mod(iyear-1,100)+1   ! two digits
 icen   = (iyear-1)/100+1
!
!KPDS     - array containing pds elements. 
!
 kpds(1) =   7             ! ID of center
 kpds(2) =  44             ! generating process ID number 
 if (imx == 4320) then
    kpds(3) = 173          ! grid definition for 1/12 deg. grid
 elseif (imx == 720) then
    kpds(3) = 235          ! grid definition for 1/2 deg. grid
 else
    kpds(3) = 255          ! catch-all grid definition
 endif

 if ( maskflag == 0 ) then
    kpds(4) = 128    ! gds/bms flag (right adj copy of octet 8)
 else
    kpds(4) = 192    ! turn on 7th bit from end for bit map
 endif

 kpds(5) =  11    ! indicator of temperature
 kpds(6) =   1    ! type of level
 kpds(7) =   0    ! height/pressure , etc of level
 kpds(8) = iyr    ! year of century
 kpds(9) = imon   ! month of year
 kpds(10)= iday   ! day of month
 kpds(11)= ihr    ! hour of day
 kpds(12)=   0    ! minute of hour
 kpds(13)=   2    ! indicator of forecast time unit
 kpds(14)=   0    ! time range 1
 kpds(15)=   0    ! time range 2
 kpds(16)=   1    ! time range flag
 kpds(17)=   0    ! number included in average
 kpds(18)=   1    ! version nr of grib specification
 kpds(19)=   3    ! version nr of parameter table
 kpds(20)=   0    ! nr missing from average/accumulation
 kpds(21)= icen   ! century of reference time of data
 kpds(22)=   2    ! units decimal scale factor
 kpds(23)=   4    ! subcenter number
 kpds(24)=   0    ! pds byte 29, for nmc ensemble products
!                   128 IF FORECAST FIELD ERROR
!                    64 IF BIAS CORRECTED FCST FIELD
!                    32 IF SMOOTHED FIELD
!                    WARNING: CAN BE COMBINATION OF MORE THAN 1
 kpds(25)=   0    !  - PDS BYTE 30, NOT USED
!       (26-35)  - RESERVED
!       (36-N)   - CONSECUTIVE BYTES EXTRACTED FROM PROGRAM
!                  DEFINITION SECTION (PDS) OF GRIB MESSAGE
!
!KGDS     - ARRAY CONTAINING GDS ELEMENTS.
 kgds(1) =     0                         ! data representation type
 kgds(2) =   imx                         ! points on latitude circle
 kgds(3) =   imy                         ! points on longitude meridian
 loninc  = nint(360000./real(imx))       ! lon increment
 latinc  = nint(180000./real(imy))       ! lat increment
!
! the latitide from S to N for 1/4 and 1 degree
! the latitide from N to S for others
!
 if ( imx == 1440 ) then
    kgds(4) = -90000+nint(90000./real(imy)) ! latitude of origin
    kgds(7) =  90000-nint(90000./real(imy)) ! latitude of extreme point
    kgds(11)=    64                         ! scanning mode flag (right adj copy of octet 28)
 else
    kgds(4) =  90000-nint(90000./real(imy)) ! latitude of origin
    kgds(7) = -90000+nint(90000./real(imy)) ! latitude of extreme point
    kgds(11)=     0                         ! scanning mode flag (right adj copy of octet 28)
 endif
 kgds(5) = nint(180000./real(imx))          ! longitude of origin
 kgds(6) =   128                            ! resolution flag (right adj copy of octet 17)
 kgds(8) = nint(-180000./real(imx))         ! longitude of extreme point
 kgds(9) = latinc                           ! di Latitudinal direction of increment
 kgds(10)= loninc                           ! dj longitudinal direction increment
 kgds(19)=     0                            ! number of vertical coordinate parameters
 kgds(20)=   255                            ! octet number of the list of vertical coordinate parameters
                                            !    OR
                                            ! octet number of the list of numbers of points in each row
                                            !    OR
                                            !  255 if neither are present
 kgds(21)=     0                            ! for grids with pl, number of points in grid
 kgds(22)=     0                            ! number of words in each row

!print*,'kpds:',kpds

 do i=12,18
   kgds(i) = 0    ! 12-18 not used
 enddo
!print*,'kgds:',kgds

!
! save data from S to N for 1/4 and 1 degree resolutions
! save data from N to S for others
!
 if ( imx == 1440 ) then
    do j = 1, imy
       do i = 1, imx
          out(i,j)=anal(i,j)
          out_mask(i,j)=mask(i,j)
       enddo
       if ( maskflag /= 0 ) then
          do i = 1, imx
             lb(i,j) = mask(i,j) == 0
          enddo
       endif
    enddo
 else
    do j = 1, imy
       jj = imy+1-j
       do i = 1, imx
          out(i,jj)=anal(i,j)
          out_mask(i,jj)=mask(i,j)
       enddo
       if ( maskflag /= 0 ) then
          do i = 1, imx
             lb(i,jj) = mask(i,j) == 0
          enddo
       endif
    enddo
 endif

 call baopen(lun_grb,fname,iret)
 call putgb(lun_grb,kf,kpds,kgds,lb,out,iret)
!
! store mask (0 or 1) to grib file for 1/4 degree OISST-like 
!
 if (imx == 1440) then
    kpds(5) =  81                       ! indicator of mask
    kpds(22)=   0                       ! units decimal scale factor
    call putgb(lun_grb,kf,kpds,kgds,lb,out_mask,iret)
    write(*,*) 'add mask to grib file'
 endif
 print*,'iret from putgb to fname for anal is ',iret,fname
 call baclose(lun_grb,iret)

 end subroutine togrib
