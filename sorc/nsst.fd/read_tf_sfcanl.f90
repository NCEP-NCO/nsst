subroutine read_tf_sfcanl_nems(filename,tf,slmsk,xlats,xlons,nlat,nlon)
!$$$  main program documentation block
!
! subroutine:  read_tf_sfcanl_nems
!
! prgmmr: Xu Li         org: noaa               date: 2019-03-13
!
! abstract:   read tf & mask from GFS FV3 SFC anal file (nemsio)
!
! program history log:
!
! attributes:
!   language: f95
!
!$$$

  use nemsio_module, only: nemsio_init,nemsio_open,nemsio_close
  use nemsio_module, only: nemsio_gfile,nemsio_getfilehead,nemsio_readrecv

  implicit none

  character(len=6), intent(in) :: filename
  integer, intent(in) :: nlon,nlat
  real,    dimension(nlon,nlat), intent(out) :: tf
  integer, dimension(nlon,nlat), intent(out) :: slmsk
  real,    dimension(nlat), intent(out) :: xlats
  real,    dimension(nlon), intent(out) :: xlons
! local 
  real(4), dimension(nlon*nlat) :: rwork1d
  real, dimension(nlon,nlat) :: rwork2d

  logical :: nemsio
  integer :: iret,nrec
  integer :: lonb, latb, i, j
  integer, dimension(7):: idate

  type(nemsio_gfile) :: gfile

  nemsio=.false.

  call nemsio_init(iret)

  call nemsio_open(gfile,trim(filename),'READ',iret)
  if (iret == 0 ) then
     nemsio = .true.
  else
     write(6,*)'***ERROR*** ',trim(filename),' contains unrecognized format.  ABORT'
  endif

  if (.not.nemsio ) then
     write(6,*)'invalid sfc format'
     stop
  endif

  call nemsio_getfilehead(gfile, nrec=nrec, idate=idate, dimx=lonb, dimy=latb, iret=iret)

  if ( lonb /= nlon .or. latb /= nlat ) then
     write(*,*) ' inconsistent dimensions, nlon, nlat = ',nlon, nlat,' vs lonb, latb = ',lonb, latb
  endif
        
  rwork1d=0.0
  rwork2d=0.0

  call nemsio_getfilehead(gfile, lon=rwork1d, iret=iret)
  rwork2d=reshape(rwork1d,(/size(rwork2d,1),size(rwork2d,2)/))
  xlons(:) = rwork2d(:,1)

  call nemsio_getfilehead(gfile, lat=rwork1d, iret=iret)
  rwork2d=reshape(rwork1d,(/size(rwork2d,1),size(rwork2d,2)/))
  do j = 1, nlat
     do i = 1, nlon
        xlats(j) = rwork2d(1,nlat+1-j)
     enddo
  enddo

  call nemsio_readrecv(gfile,'land','sfc',1,rwork1d,iret)
  rwork2d=reshape(rwork1d,(/size(slmsk,1),size(slmsk,2)/))
  do j = 1, nlat
     do i = 1, nlon
        slmsk(i,j) = int(rwork2d(i,nlat+1-j))
     enddo
  enddo

  call nemsio_readrecv(gfile,'tref', 'sfc',1,rwork1d,iret)
  rwork2d=reshape(rwork1d,(/size(tf,1),size(tf,2)/))
  do j = 1, nlat
     do i = 1, nlon
        tf(i,j) = rwork2d(i,nlat+1-j)
     enddo
  enddo

  call nemsio_close(gfile, iret)

end subroutine read_tf_sfcanl_nems

subroutine get_sfcanl_dim_nems(filename,nlat,nlon)
!
! subroutine:  sfcanl_dim_get
!
! prgmmr: Xu Li         org: noaa               date: 2019-03-13
!
! abstract:   get dimensions of nemsio sfcanl 
!
! program history log:
!
! attributes:
!   language: f95
!
!$$$

  use nemsio_module, only: nemsio_init,nemsio_open,nemsio_close
  use nemsio_module, only: nemsio_gfile,nemsio_getfilehead

  implicit none

  character(len=6), intent(in) :: filename
  integer, intent(out) :: nlat,nlon

  integer :: iret,nrec
  integer, dimension(7):: idate

  logical :: nemsio
  type(nemsio_gfile) :: gfile

  nemsio=.false.

  write(*,*) 'filename : ',trim(filename)
  call nemsio_init(iret)

  call nemsio_open(gfile,trim(filename),'READ',iret)
  if (iret == 0 ) then
     nemsio = .true.
  else
     write(6,*)'***ERROR*** ',trim(filename),' contains unrecognized format.  ABORT'
  endif

  if (.not.nemsio ) then
     write(6,*)'invalid sfc format'
     stop
  endif

  call nemsio_getfilehead(gfile, nrec=nrec, idate=idate, dimx=nlon, dimy=nlat, iret=iret)

  write(*,*) 'get_sfcanl_dim_nemsio : nlat, nlon :',nlat, nlon 

  call nemsio_close(gfile, iret)

 end subroutine get_sfcanl_dim_nems

subroutine read_tf_sfcanl_nc(filename,tf,mask,xlats,xlons,nlat,nlon)

  use netcdf
  implicit none

  ! This is the name of the data file we will read.
  character (len=6),             intent(in)  :: filename
  integer,                       intent(in)  :: nlat,nlon
  real,    dimension(nlat),      intent(out) :: xlats
  real,    dimension(nlon),      intent(out) :: xlons
  real,    dimension(nlon,nlat), intent(out) :: tf
  integer, dimension(nlon,nlat), intent(out) :: mask
! Local variables
  integer :: ncid

  integer, dimension(nlon,nlat) :: mask0
  real,    dimension(nlon,nlat) :: tf0
  real,    dimension(nlat) :: xlats0

  integer, parameter :: ndims = 3
  character (len = *), parameter :: lat_name = "grid_yt"
  character (len = *), parameter :: lon_name = "grid_xt"
  character (len = *), parameter :: mask_name="land"
  character (len = *), parameter :: tf_name="tref"
! mask = 0 (water), 1 (land), 2 (sea ice)
  integer :: lon_varid,lat_varid,tf_varid,mask_varid


  ! The start and count arrays will tell the netCDF library where to read our
  ! data.
  integer, dimension(ndims) :: start, count

  character (len = *), parameter :: units = "units"
  character (len = *), parameter :: tf_units = "kelvin", mask_units = "none"
  character (len = *), parameter :: time_units = "second"
  character (len = *), parameter :: lat_units = "degrees_north"
  character (len = *), parameter :: lon_units = "degrees_east"

! Loop indices
  integer :: i,j
  write(*,*) 'process netcdf sfc file : ',filename

! Open the file.
  call nc_check( nf90_open(filename, nf90_nowrite, ncid) )

! Get the varids of the latitude and longitude coordinate variables.
  call nc_check( nf90_inq_varid(ncid, "grid_yt", lat_varid) )
  call nc_check( nf90_inq_varid(ncid, "grid_xt", lon_varid) )

! Read the latitude and longitude data.
  call nc_check( nf90_get_var(ncid, lat_varid, xlats0) )
  call nc_check( nf90_get_var(ncid, lon_varid, xlons) )

  do j = 1, nlat
     xlats(j) = xlats0(nlat+1-j)
  enddo

! Get the varids of the tf & mask netCDF variables.
  call nc_check( nf90_inq_varid(ncid, "land",   mask_varid) )
  call nc_check( nf90_inq_varid(ncid, "tref",   tf_varid) )

! Read 1 record of nlat*nlon values, starting at the beginning
! of the record (the (1, 1, rec) element in the netCDF file).
  start = (/ 1, 1, 1 /)
  count = (/ nlon, nlat, 1 /)

! Read the tf & mask data from the file, one record at a time.
  call nc_check( nf90_get_var(ncid, mask_varid,   mask0,start, count) )
  call nc_check( nf90_get_var(ncid, tf_varid,     tf0,  start, count) )


  do j = 1, nlat
     do i = 1, nlon
        mask(i,j) = mask0(i,nlat+1-j)
        tf(i,j)   = tf0(i,nlat+1-j)
     enddo
  enddo

! Close the file. This frees up any internal netCDF resources
! associated with the file.
  call nc_check( nf90_close(ncid) )

! If we got this far, everything worked as expected. Yipee!
  print *,"*** SUCCESS reading file ", filename, "!"

end subroutine read_tf_sfcanl_nc

subroutine get_sfcanl_dim_nc(filename,nlat,nlon)
! abstract: get dimensions of sal array
  use netcdf

  character (len=*), intent(in)  :: filename
  integer,           intent(out) :: nlat,nlon
! Local variables
  character (len = *), parameter :: lat_name = "grid_yt"
  character (len = *), parameter :: lon_name = "grid_xt"
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

end subroutine get_sfcanl_dim_nc

