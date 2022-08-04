subroutine get_seaice(xlats_ij,xlons_ij,ny,nx,sice)
!
! Abstract: get sea ice analysis at atime and target resolution (nx,ny)
!
!
! input & output
!
 real,    dimension(nx*ny), intent(in)  :: xlats_ij   ! latitudes of target grids (nx*ny)
 real,    dimension(nx*ny), intent(in)  :: xlons_ij   ! latitudes of target grids (nx*ny)
 real,    dimension(nx,ny), intent(out) :: sice       ! sea ice atime (nx,ny)
 integer, intent(in) :: nx,ny
! local declare
 integer, dimension(nx*ny)              :: mask_ij    ! mask at target grids (0 = water, 1 = land, 2 = sea ice) (nx*ny)
 real,    dimension(nx*ny)              :: sice_ij    ! sea ice  analysis at target grids (nx*ny)

 real,    allocatable, dimension(:,:)   :: sice0      ! sea ice analysis
 integer, allocatable, dimension(:,:)   :: smask      ! surface mask of sea ice
 real,    allocatable, dimension(:)     :: sxlats     ! latitudes of sea ice analysis
 real,    allocatable, dimension(:)     :: sxlons     ! longitudes of sea ice analysis

 integer :: nxs,nys,mon1,mon2,sfcflag,i,j
 character (len=6), parameter :: fin_iceanl='iceanl'  ! sea ice analysis file name
!
! get the dimensions of the sea ice analysis & allocate the related arrays
!
 call get_dim_grb(fin_iceanl,nys,nxs)
 allocate( sice0(nxs,nys),smask(nxs,nys),sxlats(nys),sxlons(nxs) )
!
! read sea ice analysis
!
  call read_seaice_grb(fin_iceanl,nys,nxs,sxlats,sxlons,sice0) 
!
! Handle the filler (1.57) over land and others
!
 do j = 1, ny
    do i = 1, nx
       if ( sice0(i,j) > 1.0 .or. sice0(i,j) < 0.0 ) sice0(i,j) = 0.0
    enddo
 enddo
!
! get sice (nx by ny lat/lon)
!
 if ( nx == nxs .and. ny == nys ) then
    sice(:,:) = sice0(:,:)
    write(*,'(a,4I8,2F9.3)') 'get_seaice: same dimensions,nx,ny,nxs,nys : ',nx,ny,nxs,nys,minval(sice),maxval(sice)
 else
    sfcflag=0
    smask=0
    mask_ij=0
!   write(*,'(a,4I8)') 'different dimensions,nx,ny,nxs,nys : ',nx,ny,nxs,nys
    call lalo_to_tile(sice0,   smask,   sxlats,  sxlons,  nys, nxs, &
                      sice_ij, mask_ij, xlats_ij,xlons_ij,ny,  nx, &
                      sfcflag,0.0,0,0.0)
!   write(*,'(a,2F9.3)') 'done with lalo_to_tile for sice',minval(sice_ij),maxval(sice_ij)
    sice(:,:) = reshape (sice_ij, (/nx,ny/) )
 endif
end subroutine get_seaice

subroutine get_dim_grb(file_ice,mlat,mlon)
!                .      .    .                                       .
! abstract:   get_grb_dim :  get dimensions
!   prgmmr: xu li            org: np23                date: 2019-03-13
!
!
! program history log:
!
!   input argument list:
!     file - file name of GRIB SST file
! output
!     mlat,mlon
!     call subs: getgbh
!
!
!$$$
  implicit none

! Declare passed variables and arrays
  character (len=*),  intent(in)  :: file_ice
  integer,            intent(out) :: mlat,mlon
! Declare local parameters
  integer,parameter:: lu = 22   ! FORTRAN unit number of GRIB SST file

  integer :: iret,i,j,kf,kg,k
  integer, dimension(22):: jgds,kgds
  integer, dimension(25):: jpds,kpds

!************+******************************************************************************
!
! Open analysis file (GRIB)
  write(*,*) 'get_dim_grb, file_ice : ',file_ice
  call baopenr(lu,trim(file_ice),iret)
  if (iret /= 0 ) then
     write(6,*)'get_dim_grb:  ***ERROR*** opening sea ice file',file_ice
     stop
  endif

! Define variables for read
  j=-1
  jpds=-1
  jgds=-1
  jpds(5)=91          ! sea ice variable
  jpds(6)=102         ! surface
  call getgbh(lu,0,j,jpds,jgds,kg,kf,k,kpds,kgds,iret)

  mlat = kgds(3)   ! number points on longitude circle (360)
  mlon = kgds(2)   ! number points on latitude circle (720)

  write(*,*) 'sea ice nys, nxs : ',mlat, mlon

  call baclose(lu,iret)
  if (iret /= 0 ) then
     write(6,*)'get_dim_grb:  ***ERROR*** close sea ice file'
     stop
  endif

end subroutine get_dim_grb

subroutine read_seaice_grb(file_sice,mlat,mlon,xlats,xlons,sice)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    rdgrbsice                   read grib1 SST analysis
!   prgmmr: xu li            org: np23                date: 2003-04-04
!
! abstract: read sea ice analysis (GRIB format) and save it as expanded and transposed array
!
!     Subroutine rdgrbsice must be compiled with the NCEP W3 library
!     and the BACIO library.
!
!
! program history log:
!   2003-04-04  xu li, bert katz
!   2005-04-18  treadon - fill southern and northern rows of sice grid
!                         with mean of adjacent row.  This treatment is
!                         consistent with rdgesfc.f90 and rdgesig.f90
!
!   input argument list:
!     file_sice - file name of GRIB sea ice file
!     mlat      - dimension in y-direction
!     mlon      - dimension in x-direction
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
  character(6),               intent(in)  :: file_sice
  integer,                    intent(in)  :: mlat,mlon
  real, dimension(mlat),      intent(out) :: xlats
  real, dimension(mlon),      intent(out) :: xlons
  real, dimension(mlon,mlat), intent(out) :: sice

! Declare local parameters
  integer,parameter:: lu_sice = 22   ! FORTRAN unit number of GRIB SST file
  real, parameter :: deg2rad = 3.141593/180.0

! Declare local variables and arrays
  logical(1), allocatable, dimension(:) ::  lb

  integer :: nlat_sice,nlon_sice,ilon,jlat
  integer :: iret,ni,nj
  integer :: mscan,kb1
  integer :: jincdir,i,iincdir,kb2,kb3,kf,kg,k,j,jf
  integer, dimension(22):: jgds,kgds
  integer, dimension(25):: jpds,kpds

  real :: x0,y0,dres
  real, allocatable, dimension(:) :: f
  
!************+******************************************************************************
!
! Open SST analysis file (GRIB)
  call baopenr(lu_sice,trim(file_sice),iret)
  if (iret /= 0 ) then
     write(6,*)'RDGRBICE:  ***ERROR*** opening Sea Ice file'
     stop
  endif
! Define SST variables for read
  j=-1
  jpds=-1
  jgds=-1
  jpds(5)=91          ! Sea Ice fraction (1=ice; 0=no ice)
! jpds(6)=102         ! surface
  call getgbh(lu_sice,0,j,jpds,jgds,kg,kf,k,kpds,kgds,iret)

  nlat_sice = kgds(3)   ! number points on longitude circle (360)
  nlon_sice = kgds(2)   ! number points on latitude circle (720)

  if ( mlat /= nlat_sice .or. mlon /= nlon_sice ) then
     write(*,*) 'STOP inconsistent dimensions : ',nlat_sice, nlon_sice,mlat,mlon
     stop
  endif
!
! get xlats and xlons
!
  dres = 180.0/real(mlat)
  y0 = 0.5*dres-90.0
  x0 = 0.5*dres

! Get lat_sst & lon_sst
  do jlat = 1, mlat
     xlats(jlat) = y0 + real(jlat-1)*dres
  enddo

  do ilon = 1, mlon
     xlons(ilon) = (x0 + real(ilon-1)*dres)
  enddo
! Allocate arrays
  allocate(lb(nlat_sice*nlon_sice))
  allocate(f(nlat_sice*nlon_sice))
  jf=nlat_sice*nlon_sice

! Read in the analysis
  call getgb(lu_sice,0,jf,j,jpds,jgds,kf,k,kpds,kgds,lb,f,iret)
  if (iret /= 0) then
     write(6,*)'RDGRBSICE:  ***ERROR*** reading sice analysis data record'
     deallocate(lb,f)
     stop
  endif

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
     sice(i,j)=f(k)
  end do
     
  deallocate(lb,f)
  
end subroutine read_seaice_grb
