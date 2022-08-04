subroutine get_tf_clm(mask_ij,xlats_ij,xlons_ij,ny,nx,iy,im,id,ih,dsearch,miss_fill,bmiss,tf_clm)
!
! Abstract: get SST climatology at the valid time (atime) and target resolution (nx,ny)
!
! input & output 
!
 real,    dimension(nx*ny), intent(in)  :: xlats_ij   ! latitudes of target grids (nx*ny)
 real,    dimension(nx*ny), intent(in)  :: xlons_ij   ! latitudes of target grids (nx*ny)
 integer, dimension(nx*ny), intent(in)  :: mask_ij    ! mask at target grids (0 =water, 1 = land, 2 = sea ice) (nx*ny)
 real,    dimension(nx,ny), intent(out) :: tf_clm     ! SST climatology valid at atime (nx,ny)
 real,    intent(in) :: dsearch,bmiss
 integer, intent(in) :: iy,im,id,ih,nx,ny,miss_fill
! local declare
 real,    allocatable, dimension(:,:)   :: tf_clm0    ! SST climatology at the valid time
 integer, allocatable, dimension(:,:)   :: cmask      ! surface mask of SST climatology (set to zero = water everywhere)
 real,    allocatable, dimension(:)     :: cxlats     ! latitudes of SST Climatology
 real,    allocatable, dimension(:)     :: cxlons     ! longitudes of SST Climatology

 real,    dimension(nx*ny)  :: tf_clm_ij  ! SST climatology at target grids (nx*ny)
 integer, dimension(nx*ny)  :: cmask_ij  ! SST climatology at target grids (nx*ny)
 integer :: nxc,nyc,mon1,mon2,sfcflag,i,j
 character (len=6), parameter :: fin_tf_clm='sstclm' ! SST climatology file name
!
! get which two months used and their weights from atime
!
 call get_tim_wei(iy,im,id,ih,mon1,mon2,wei1,wei2)
!
! get the dimensions of the SST climatology & allocate the related arrays
!
 call get_tf_clm_dim(fin_tf_clm,nyc,nxc)
 allocate( tf_clm0(nxc,nyc),cmask(nxc,nyc),cxlats(nyc),cxlons(nxc) )
!
! get tf_clm at the analysis time from monthly climatology & cxlats, cxlons
!
 call get_tf_clm_ta(tf_clm0,cxlats,cxlons,nyc,nxc,mon1,mon2,wei1,wei2)
!
! get tf_clm (nx by ny lat/lon) valid at atime
!
 if ( nx == nxc .and. ny == nyc ) then
    tf_clm(:,:) = tf_clm0(:,:)
    write(*,'(a,2F9.3)') 'same dimensions, tf_clm, min, max : ',minval(tf_clm),maxval(tf_clm)
 else
    write(*,'(a,4I8)') ' get_tf_clm: different dimensions,nx,ny,nxc,nyc : ',nx,ny,nxc,nyc
    sfcflag=0
    cmask = 0
    cmask_ij = 0
    call lalo_to_tile(tf_clm0,   cmask,    cxlats,  cxlons,  nyc, nxc, &
                      tf_clm_ij, cmask_ij, xlats_ij,xlons_ij,ny,  nx,  &
                      sfcflag,dsearch,miss_fill,bmiss)
    write(*,'(a,2F9.3)') 'tf_clm0, min, max                           : ',minval(tf_clm0),maxval(tf_clm0)
    write(*,'(a,2F9.3)') 'done with lalo_to_tile for tf_clm, min, max : ',minval(tf_clm_ij),maxval(tf_clm_ij)

    tf_clm(:,:) = reshape (tf_clm_ij, (/nx,ny/) )
 endif

end subroutine get_tf_clm

subroutine read_tf_clim_grb(file_sst,sst,rlats_sst,rlons_sst,mlat_sst,mlon_sst,mon)

!                .      .    .                                       .
! abstrac:    read_tf_clim_grb :  read grib1 SST analysis
!   prgmmr: xu li            org: np23                date: 2019-03-13
!
! abstract: read SST analysis (GRIB format) and save it as expanded and transposed array
!
!     Subroutine rdgrbsst must be compiled with the NCEP W3 library
!     and the BACIO library.
!
!
! program history log:
!
!   input argument list:
!     file_sst - file name of GRIB SST file
!     mon      - month number
!     mlat_sst,mlon_sst
!   output:
!     rlats_sst
!     rlons_sst
!     sst
!
!   argument list defined by this reading:
!     sst - SST field
!     Note: (1) The data is stored from north to south originally in GRIB format,
!             but is stored from south to north with this reading routine
!     nlat_sst  - latitudinal dimension of SST
!     nlon_sst  - longitudinal dimension of SST
!     xsst0 - latitude of origin
!     ysst0 - longitude of origin
!     dres  - lat/lon increment
!
!     call subs: getgbh, getgb
!
! attributes:
!   language: f90
!
!$$$
  implicit none

! Declare passed variables and arrays
  character(*)                       , intent(in   ) :: file_sst
  integer                            , intent(in   ) :: mon,mlat_sst,mlon_sst
  real, dimension(mlat_sst)          , intent(  out) :: rlats_sst
  real, dimension(mlon_sst)          , intent(  out) :: rlons_sst
  real, dimension(mlon_sst,mlat_sst) , intent(  out) :: sst

! Declare local parameters
  integer,parameter:: lu_sst = 21   ! FORTRAN unit number of GRIB SST file
  real, parameter :: deg2rad = 3.141593/180.0

! Declare local variables and arrays
  logical(1), allocatable, dimension(:)    ::  lb

  integer :: nlat_sst,nlon_sst
  integer :: iret,ni,nj
  integer :: mscan,kb1
  integer :: jincdir,i,iincdir,kb2,kb3,kf,kg,k,j,jf
  integer, dimension(22):: jgds,kgds
  integer, dimension(25):: jpds,kpds

  real :: xsst0,ysst0,dres
  real, allocatable, dimension(:) :: f
  
!************+******************************************************************************
!
! Open SST analysis file (GRIB)
  write(*,*) ' sstclm : ',file_sst
  call baopenr(lu_sst,trim(file_sst),iret)
  if (iret /= 0 ) then
     write(6,*)'read_tf_clm_grb:  ***ERROR*** opening SST file'
     stop
  endif

! Define SST variables for read
  j=-1
  jpds=-1
  jgds=-1
  jpds(5)=11        ! SST variable
  jpds(6)=1         ! surface
  jpds(9) = mon
  call getgbh(lu_sst,0,j,jpds,jgds,kg,kf,k,kpds,kgds,iret)

  nlat_sst = kgds(3)   ! number points on longitude circle (360)
  nlon_sst = kgds(2)   ! number points on latitude circle (720)

! write(*,*) 'nlat_sst, nlon_sst, mon : ',nlat_sst, nlon_sst, mon
! write(*,*) 'mlat_sst, mlon_sst : ',mlat_sst, mlon_sst

! Allocate arrays
  allocate(lb(nlat_sst*nlon_sst))
  allocate(f(nlat_sst*nlon_sst))
  jf=nlat_sst*nlon_sst

! Read in the analysis
  call getgb(lu_sst,0,jf,j,jpds,jgds,kf,k,kpds,kgds,lb,f,iret)
  if (iret /= 0) then
     write(6,*)'read_tf_clm_grb:  ***ERROR*** reading sst analysis data record'
     deallocate(lb,f)
     stop
  endif

  if ( (nlat_sst /= mlat_sst) .or. (nlon_sst /= mlon_sst) ) then
     write(6,*)'read_rtg_org:  inconsistent dimensions.  mlat_sst,mlon_sst=',&
          mlat_sst,mlon_sst,' -versus- nlat_sst,nlon_sst=',nlat_sst,nlon_sst
     deallocate(lb,f)
     stop
  endif

!
! get xlats and xlons
!
  dres = 180.0/real(nlat_sst)
  ysst0 = 0.5*dres-90.0
  xsst0 = 0.5*dres

! Get lat_sst & lon_sst
  do j = 1, nlat_sst
     rlats_sst(j) = ysst0 + real(j-1)*dres
  enddo

  do i = 1, nlon_sst 
     rlons_sst(i) = (xsst0 + real(i-1)*dres)
  enddo

! Load dimensions and grid specs.  Check for unusual values
  ni=kgds(2)                ! 720 for 0.5 x 0.5
  nj=kgds(3)                ! 360 for 0.5 x 0.5 resolution

  mscan=kgds(11)
  kb1=ibits(mscan,7,1)   ! i scan direction
  kb2=ibits(mscan,6,1)   ! j scan direction
  kb3=ibits(mscan,5,1)   ! (i,j) or (j,i)

! Get i and j scanning directions from kb1 and kb2.
! 0 yields +1, 1 yields -1. +1 is west to east, -1 is east to west.
  iincdir = 1-kb1*2

! 0 yields -1, 1 yields +1. +1 is south to north, -1 is north to south.
  jincdir = kb2*2 - 1

! write(*,*) 'read_tf_clim_grb,iincdir,jincdir : ',iincdir,jincdir
  do k=1,kf

!    kb3 from scan mode indicates if i points are consecutive
!    or if j points are consecutive
     if(kb3==0)then     !  (i,j)
        i=(ni+1)*kb1+(mod(k-1,ni)+1)*iincdir
        j=(nj+1)*(1-kb2)+jincdir*((k-1)/ni+1)
     else                !  (j,i)
        j=(nj+1)*(1-kb2)+(mod(k-1,nj)+1)*jincdir
        i=(ni+1)*kb1+iincdir*((k-1)/nj+1)
     endif
     sst (i,j) = f(k)
  end do
     
  deallocate(lb,f)

  call baclose(lu_sst,iret)
  if (iret /= 0 ) then
     write(6,*)'read_tf_clm_grb:  ***ERROR*** close SST file'
     stop
  endif
  
end subroutine read_tf_clim_grb


subroutine read_msk_clim_grb(file_sst,slmsk,mlat_sst,mlon_sst)

!                .      .    .                                       .
! abstrac:    read_msk_clim_grb :  read grib1 bitmap mask
!   prgmmr: xu li            org: np23                date: 2019-03-13
!
!
! program history log:
!
!   input argument list:
!     file_sst - file name of GRIB SST file
!     mlat_sst,mlon_sst
!   output:
!     slmsk
!
!     call subs: getgbh, getgb
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  implicit none

! Declare passed variables and arrays
  character(*)                                   , intent(in   ) :: file_sst
  integer                                        , intent(in   ) :: mlat_sst,mlon_sst
  integer,          dimension(mlon_sst,mlat_sst) , intent(  out) :: slmsk

! Declare local parameters
  integer,parameter:: lu_sst = 21   ! FORTRAN unit number of GRIB SST file

! Declare local variables and arrays
  real, allocatable, dimension(:) :: f
  logical(1), allocatable, dimension(:)    ::  lb
  logical(1), dimension(mlon_sst,mlat_sst) ::  lbms

  integer :: nlat_sst,nlon_sst
  integer :: iret,ni,nj
  integer :: mscan,kb1
  integer :: jincdir,i,iincdir,kb2,kb3,kf,kg,k,j,jf
  integer, dimension(22):: jgds,kgds
  integer, dimension(25):: jpds,kpds

!************+******************************************************************************
!
! Open SST analysis file (GRIB)
  call baopenr(lu_sst,trim(file_sst),iret)
  if (iret /= 0 ) then
     write(6,*)'read_tf_clm_grb:  ***ERROR*** opening SST file'
     stop
  endif

! Define SST variables for read
  j=-1
  jpds=-1
  jgds=-1
  jpds(5)=11        ! SST variable
  jpds(6)=1         ! surface
  jpds(9) = 1       ! use month 1 since mask is not month dependent in Clim file
  call getgbh(lu_sst,0,j,jpds,jgds,kg,kf,k,kpds,kgds,iret)

  nlat_sst = kgds(3)   ! number points on longitude circle (360)
  nlon_sst = kgds(2)   ! number points on latitude circle (720)

! Allocate arrays
  allocate(lb(nlat_sst*nlon_sst))
  allocate(f(nlat_sst*nlon_sst))
  jf=nlat_sst*nlon_sst

! Read in the analysis
  call getgb(lu_sst,0,jf,j,jpds,jgds,kf,k,kpds,kgds,lb,f,iret)
  if (iret /= 0) then
     write(6,*)'read_msk_clm_grb:  ***ERROR*** reading sst analysis data record'
     deallocate(lb,f)
     stop
  endif

  if ( (nlat_sst > mlat_sst) .or. (nlon_sst > mlon_sst) ) then
     write(6,*)'read_rtg_org:  inconsistent dimensions.  mlat_sst,mlon_sst=',&
          mlat_sst,mlon_sst,' -versus- nlat_sst,nlon_sst=',nlat_sst,nlon_sst
     deallocate(lb)
     stop
  endif

! Load dimensions and grid specs.  Check for unusual values
  ni=kgds(2)               
  nj=kgds(3)              

  mscan=kgds(11)
  kb1=ibits(mscan,7,1)   ! i scan direction
  kb2=ibits(mscan,6,1)   ! j scan direction
  kb3=ibits(mscan,5,1)   ! (i,j) or (j,i)

! Get i and j scanning directions from kb1 and kb2.
! 0 yields +1, 1 yields -1. +1 is west to east, -1 is east to west.
  iincdir = 1-kb1*2

! 0 yields -1, 1 yields +1. +1 is south to north, -1 is north to south.
  jincdir = kb2*2 - 1
  do k=1,kf

!    kb3 from scan mode indicates if i points are consecutive
!    or if j points are consecutive
     if(kb3==0)then     !  (i,j)
        i=(ni+1)*kb1+(mod(k-1,ni)+1)*iincdir
        j=(nj+1)*(1-kb2)+jincdir*((k-1)/ni+1)
     else                !  (j,i)
        j=(nj+1)*(1-kb2)+(mod(k-1,nj)+1)*jincdir
        i=(ni+1)*kb1+iincdir*((k-1)/nj+1)
     endif
     lbms(i,j) = lb(k)
  end do
     
  slmsk=1
  do j = 1, mlat_sst
     do i = 1, mlon_sst
        if ( lbms(i,j) ) then
           slmsk(i,j) = 0
        endif
     enddo
  enddo

  deallocate(lb)

  call baclose(lu_sst,iret)
  if (iret /= 0 ) then
     write(6,*)'read_msk_clm_grb:  ***ERROR*** close SST Clim file'
     stop
  endif
  
end subroutine read_msk_clim_grb

subroutine get_tf_clm_dim(file_sst,mlat_sst,mlon_sst)
!                .      .    .                                       .
! abstract:   get_tf_clm_dim :  get dimension of RTG SST climatology
!   prgmmr: xu li            org: np23                date: 2019-03-13
!
!
! program history log:
!
!   input argument list:
!     file_sst - file name of GRIB SST file
! output
!     mlat_sst,mlon_sst
!
!     call subs: getgbh
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  implicit none

! Declare passed variables and arrays
  character(*)                                   , intent(in ) :: file_sst
  integer                                        , intent(out) :: mlat_sst,mlon_sst

! Declare local parameters
  integer,parameter:: lu_sst = 21   ! FORTRAN unit number of GRIB SST file

  integer :: iret
  integer :: mscan,kb1
  integer :: kf,kg,k,j
  integer, dimension(22):: jgds,kgds
  integer, dimension(25):: jpds,kpds

!************+******************************************************************************
!
! Open SST analysis file (GRIB)
  call baopenr(lu_sst,trim(file_sst),iret)
  if (iret /= 0 ) then
     write(6,*)'get_tf_clm_dim:  ***ERROR*** opening SST file'
     stop
  endif

! Define SST variables for read
  j=-1
  jpds=-1
  jgds=-1
  jpds(5)=11        ! SST variable
  jpds(6)=1         ! surface
  jpds(9) = 1
  call getgbh(lu_sst,0,j,jpds,jgds,kg,kf,k,kpds,kgds,iret)

  mlat_sst = kgds(3)   ! number points on longitude circle (360)
  mlon_sst = kgds(2)   ! number points on latitude circle (720)

  write(*,*) 'mlat_sst, mlon_sst : ',mlat_sst, mlon_sst

  call baclose(lu_sst,iret)
  if (iret /= 0 ) then
     write(6,*)'get_tf_clm_dim:  ***ERROR*** close SST file'
     stop
  endif
end subroutine get_tf_clm_dim

subroutine get_tf_clm_ta(tf_clm_ta,xlats,xlons,nlat,nlon,mon1,mon2,wei1,wei2)
!$$$
! Abstract:  get Tf/SST climatology at analysis time
! Created by Xu Li, March, 2019

 implicit none

! Input
 integer, intent(in) :: nlat,nlon,mon1,mon2
 real,    intent(in) :: wei1,wei2
! Output
 real, dimension(nlon,nlat), intent(out)   :: tf_clm_ta
 real, dimension(nlat),      intent(out)   :: xlats
 real, dimension(nlon),      intent(out)   :: xlons

!input/output data file names
 character (len=6),  parameter :: fin_tf_clm='sstclm'

! Local declare
 real, dimension(nlon,nlat) :: tf_clm1,tf_clm2

!
! read in RTG SST climatology without bitmap (surface mask) for mon1 and mon2
!
  call read_tf_clim_grb(trim(fin_tf_clm),tf_clm1,xlats,xlons,nlat,nlon,mon1)
  call read_tf_clim_grb(trim(fin_tf_clm),tf_clm2,xlats,xlons,nlat,nlon,mon2)
!
!  tf_clim at the analysis time
!
   tf_clm_ta(:,:) = wei1*tf_clm1(:,:)+wei2*tf_clm2(:,:)

 end subroutine get_tf_clm_ta

