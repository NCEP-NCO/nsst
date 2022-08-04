subroutine get_sal_clm(xlats_ij,xlons_ij,ny,nx,iy,im,id,ih,sal_clm)
!
! abstract: get salinity climatology at the valid time (atime) and target resolution (nx,ny)
!
! input & output 
!
 implicit none

 real,    dimension(nx*ny), intent(in)  :: xlats_ij   ! latitudes of target grids (nx*ny)
 real,    dimension(nx*ny), intent(in)  :: xlons_ij   ! latitudes of target grids (nx*ny)
 real,    dimension(nx,ny), intent(out) :: sal_clm    ! salinity climatology valid at atime (nx,ny)
 integer, intent(in) :: iy,im,id,ih,nx,ny
! local declare
 real,    allocatable, dimension(:,:)   :: sal_clm0   ! salinity climatology at the valid time
 real,    allocatable, dimension(:)     :: cxlats     ! latitudes of sst climatology
 real,    allocatable, dimension(:)     :: cxlons     ! longitudes of sst climatology

 real,    dimension(nx*ny)  :: sal_clm_ij  ! salinity climatology at target grids (nx*ny)
 real :: wei1,wei2
 integer :: nxc,nyc,mon1,mon2,i,j
 character (len=6), parameter :: fin_sal_clm='salclm' ! salinity climatology file name
!
! get which two months used and their weights from atime
!
 call get_tim_wei(iy,im,id,ih,mon1,mon2,wei1,wei2)
!
! get the dimensions of the sst climatology & allocate the related arrays
!
 call get_dim_nc(fin_sal_clm,nyc,nxc)
 allocate( sal_clm0(nxc,nyc),cxlats(nyc),cxlons(nxc) )
!
! get sal_clm at the analysis time from monthly climatology & cxlats, cxlons
!
 call get_sal_clm_ta(sal_clm0,cxlats,cxlons,nyc,nxc,mon1,mon2,wei1,wei2)
!
! get sal_clm (nx by ny lat/lon) valid at atime
!
 if ( nx == nxc .and. ny == nyc ) then
    sal_clm(:,:) = sal_clm0(:,:)
!   write(*,'(a,2f9.3)') 'same dimensions, sal_clm, min : ',minval(sal_clm),maxval(sal_clm)
 else
    write(*,'(a,4i8)') 'different dimensions,nx,ny,nxc,nyc : ',nx,ny,nxc,nyc
    call intp_tile(sal_clm0,  cxlats,  cxlons,  nyc, nxc, &
                   sal_clm_ij,xlats_ij,xlons_ij,ny,  nx)
!   write(*,'(a,2f9.3)') 'sal_clm0, min, max                        : ',minval(sal_clm0),maxval(sal_clm0)
!   write(*,'(a,2f9.3)') 'done with intp_tile for sal_clm, min, max : ',minval(sal_clm_ij),maxval(sal_clm_ij)

    sal_clm(:,:) = reshape (sal_clm_ij, (/nx,ny/) )
 endif

end subroutine get_sal_clm

subroutine get_sal_clm_ta(sal_clm_ta,xlats,xlons,nlat,nlon,mon1,mon2,wei1,wei2)
!$$$
! abstract:  get sal climatology at analysis time
! created by xu li, march, 2019

 implicit none

! input
 integer, intent(in) :: nlat,nlon,mon1,mon2
 real,    intent(in) :: wei1,wei2
! output
 real, dimension(nlon,nlat), intent(out)   :: sal_clm_ta
 real, dimension(nlat),      intent(out)   :: xlats
 real, dimension(nlon),      intent(out)   :: xlons

!input/output data file names
 character (len=6),  parameter :: fin_sal_clm='salclm'

! local declare
 real, dimension(nlon,nlat) :: sal_clm1,sal_clm2

!
! read in rtg sst climatology without bitmap (surface mask) for mon1 and mon2
!
  call read_salclm_gfs_nc(trim(fin_sal_clm),sal_clm1,xlats,xlons,nlat,nlon,mon1)
  call read_salclm_gfs_nc(trim(fin_sal_clm),sal_clm2,xlats,xlons,nlat,nlon,mon2)
!
!  sal_clim at the analysis time
!
   sal_clm_ta(:,:) = wei1*sal_clm1(:,:)+wei2*sal_clm2(:,:)
   write(*,'(a,2f9.3)') 'sal_clm_ta, min, max : ',minval(sal_clm_ta),maxval(sal_clm_ta)
 end subroutine get_sal_clm_ta

subroutine read_salclm_gfs_nc(filename,sal,xlats,xlons,nlat,nlon,itime)
! abstract: read woa05 salinity monthly climatology  (netcdf)
  use netcdf
  implicit none
  
  ! This is the name of the data file we will read.
  character (len=*),             intent(in)  :: filename
  integer,                       intent(in)  :: nlat,nlon
  integer,                       intent(in)  :: itime
  real,    dimension(nlat),      intent(out) :: xlats
  real,    dimension(nlon),      intent(out) :: xlons
  real,    dimension(nlon,nlat), intent(out) :: sal
! Local variables
  integer :: ncid,ntime

  integer, parameter :: ndims = 3
  character (len = *), parameter :: lat_name = "latitude"
  character (len = *), parameter :: lon_name = "longitude"
  character (len = *), parameter :: t_name = "time"
  character (len = *), parameter :: sal_name="sal"
  integer :: no_fill,fill_value
  integer :: time_varid,lon_varid, lat_varid, z_varid, sal_varid

  ! The start and count arrays will tell the netCDF library where to read our data.
  integer, dimension(ndims) :: start, count

  character (len = *), parameter :: units = "units"
  character (len = *), parameter :: sal_units = "psu" 
                                  ! PSU (Practical SalinitUnit). 1 PSU = 1g/kg
  character (len = *), parameter :: time_units = "months"
  character (len = *), parameter :: lat_units = "degrees_north"
  character (len = *), parameter :: lon_units = "degrees_east"

  integer :: missv
! Loop indices
  integer :: i,j

! Open the file. 
  call nc_check( nf90_open(filename, nf90_nowrite, ncid) )

! Get the varids of time, latitude, longitude & depth coordinate variables.
  call nc_check( nf90_inq_varid(ncid, t_name,   time_varid) )
  call nc_check( nf90_inq_varid(ncid, lat_name, lat_varid) )
  call nc_check( nf90_inq_varid(ncid, lon_name, lon_varid) )

! Read the time, latitude and longitude data.
! call nc_check( nf90_get_var(ncid, time_varid, ntime) )
  call nc_check( nf90_get_var(ncid, lat_varid,  xlats) )
  call nc_check( nf90_get_var(ncid, lon_varid,  xlons) )

! Get the varids of the sal netCDF variables.
  call nc_check( nf90_inq_varid(ncid, sal_name,sal_varid) )

! Read 1 record of nlat*nlon values, starting at the beginning 
! of the record (the (1, 1, 1, rec) element in the netCDF file).
  start = (/ 1, 1, itime /)
  count = (/ nlon, nlat, 1 /)

!  write(*,*) 'read_salclm_gfs_nc itime : ',itime
! Read the sal data from the file, one record at a time.
  call nc_check( nf90_get_var(ncid, sal_varid, sal, start, count) )

! Close the file. This frees up any internal netCDF resources
! associated with the file.
  call nc_check( nf90_close(ncid) )

! If we got this far, everything worked as expected. Yipee! 
!  print *,"*** SUCCESS reading file ", filename, "!"

end subroutine read_salclm_gfs_nc

subroutine get_dim_nc(filename,nlat,nlon)
! abstract: get dimensions of sal array
  use netcdf
  
  character (len=*), intent(in)  :: filename
  integer,           intent(out) :: nlat,nlon
! Local variables
  character (len = *), parameter :: lat_name = "latitude"
  character (len = *), parameter :: lon_name = "longitude"
  integer :: ncid
  integer :: LatDimID,LonDimID

! Open the file. 
  call nc_check( nf90_open(filename, nf90_nowrite, ncid) )

! Get dimensions 
  call nc_check( nf90_inq_dimid(ncid,lat_name,LatDimID) )
  call nc_check( nf90_inq_dimid(ncid,lon_name,LonDimID) )
  call nc_check( nf90_inquire_dimension(ncid,LatDimID,len=nlat) )
  call nc_check( nf90_inquire_dimension(ncid,LonDimID,len=nlon) )

!  write(*,'(a,1x,a6,2I8)') 'get_dim_nc, file, nlat, nlon : ',filename,nlat,nlon

! Close the file. This frees up any internal netCDF resources
! associated with the file.
  call nc_check( nf90_close(ncid) )

! If we got this far, everything worked as expected. Yipee! 
!  print *,"*** SUCCESS get dimensions from nc file ", filename, "!"

end subroutine get_dim_nc

