program nsst
!
!$$$  main program documentation block
! Abstract:  Generate GRIB SST analysis files at lat/lon grids
!           (resolution defined by namelist variables).
!
!            The sst array is from south to north n the process 
!            but it is writen from north to south in grib format
!           
! Notes of the procudure to generate RTG-like SST file from NSST Tf analysis :
!       0. Lat: South to North; Lon: from 0.5*dres to east
!       1. Read namelist variables
!       2. Get lat/lon for the target grids (ny) (nx)
!       3. Read surface mask at the target grids (nx,ny)
!       4. Get 6-hourly GFS NSST foundation temperature analysis (T1534 Gaussian grids, at prsent))
!       5. Get sea ice analysis at target grids (interpolated when necessary)
!       6. Modify the mask with sea ice info
!       7. get tf_clm at the analysis time & at target resolution from monthly climatology (nx,ny)
!       8. get sal_clm at the analysis time & at target resolution from monthly climatology (nx,ny)
!       9. interpolate tf_nst (Gaussian lat/lon) to tf (lat/lon)
!      10. set tf to be SST climatology for the missing water grids
!      11. get water temperature for grids with positive sea ice fraction 
!      12. fill the land grids temperature with known water grids temperatures
!      13. Write tf as grib format w/o bitmap
!      14. Generate 0.5 degree files while processing 1/12 degree files
!      15. Generate 1.0 degree files as gdas_oisst while processing 1/12 degree files
!     
! Created by Xu Li, Mar., 2019

!$$$

 use set_para, only: init_setup,sfcio

 implicit none

 integer, parameter :: nx_0p5=720,ny_0p5=360
 integer, parameter :: nx_1p0=360,ny_1p0=180
 integer, parameter :: sfcflag=0,miss_fill=0
 real,    parameter :: bmiss=-999.0, tzero=273.16,tfrozen=271.20
 real,    dimension(nx_0p5,ny_0p5) :: tf_0p5                  ! half degree SST aray
 integer, dimension(nx_0p5,ny_0p5) :: mask_0p5                ! half degree mask aray
 real,    dimension(nx_1p0,ny_1p0) :: tf_1p0,rmask_1p0        ! one degree SST & real mask aray
 integer, dimension(nx_1p0,ny_1p0) :: mask_1p0                ! one degree integer mask aray
!
! Arrays at Gaussian grids
!
 real,    allocatable, dimension(:,:) :: tf_nst               ! NSST Tf analysis from nemsio file (Gaussian grids)
 integer, allocatable, dimension(:,:) :: mask_nst             ! mask of NSST Tf (0 = water, others = non-water)
 real,    allocatable, dimension(:)   :: xlats_nst            ! latitudes  of NSST Tf at Gaussian grids
 real,    allocatable, dimension(:)   :: xlons_nst            ! longitudes of NSST Tf at Gaussian grids
!
! Arrays at target grids (nx,ny), (nx) or (ny)
!
 real,    allocatable, dimension(:,:) :: tf                   ! NSST Tf at target grids (nx,ny)
 real,    allocatable, dimension(:,:) :: sice                 ! Sea ice analysis at target grids (nx,ny)
 real,    allocatable, dimension(:,:) :: sal_clm              ! Salinity Climatology at target grids (nx,ny)
 real,    allocatable, dimension(:,:) :: tf_clm               ! SST Climatology at target grids and valid time (nx,ny)
 real,    allocatable, dimension(:)   :: xlats                ! latitudes of target grids (ny)
 real,    allocatable, dimension(:)   :: xlons                ! longitudes of target grids (nx)
 integer, allocatable, dimension(:,:) :: mask                 ! mask at target grids (0 = water, 1 = others) (nx,ny)
 integer, allocatable, dimension(:,:) :: masks                ! mask at target grids (0 = water, 1 = land, 2 = sea ice) (nx,ny)
!
! Arrays target grids (nx*ny)
!
 real,    allocatable, dimension(:)   :: tf_ij                ! NSST Tf at target grids (nx*ny)
 real,    allocatable, dimension(:)   :: xlats_ij             ! latitudes of target grids (nx*ny)
 real,    allocatable, dimension(:)   :: xlons_ij             ! longitudes of target grids (nx*ny)
 integer, allocatable, dimension(:)   :: mask_ij              ! mask at target grids (0 =water, 1 = land, 2 = sea ice) (nx*ny)

 real :: tfreez

!input/output data file names
 character (len=6),  parameter :: fin_tf_clm='sstclm'         ! SST climatology 
 character (len=6),  parameter :: fin_sfcanl='sfcanl'         ! GFS surface analysis file, 6-hourly
 character (len=6),  parameter :: fin_iceanl='iceanl'         ! NCEP se ice analysis file, daily
 character (len=6),  parameter :: fin_salclm='salclm'         ! Levitus Salinity Climatology (one degree)
 character (len=6),  parameter :: fin_grbmsk='grbmsk'         ! land/water mask
 character (len=10), parameter :: fin_grbmsk_0p5='grbmsk_0p5' ! land/water mask at 0.5 degree


 character (len=6), parameter :: fout_tf_grb='tf_grb'
 character (len=12), parameter :: fout_tf_grb_awips='tf_grb_awips'

 character (len=10), parameter :: fout_tf_grb_0p5='tf_grb_0p5'
 character (len=16), parameter :: fout_tf_grb_0p5_awips='tf_grb_0p5_awips'

 character (len=10), parameter :: fout_tf_grb_1p0='tf_grb_1p0'

 integer :: lensfc,mon1,mon2
 integer :: iy,im,id,ih,i,j,k,ij,ii,jj,ix,jx
 integer :: nxc,nyc,nxc_sal,nyc_sal
 integer :: nx_nst,ny_nst,nx,ny
 integer :: nfill_clm_water,nfill_clm_land,nfill_clm_ice
 real    :: wei1,wei2,dsearch,x0,y0,dres
 integer :: maskflag      ! 1 = with bitmap; 0 = w/o bitmap
 logical :: lputsi
 character (len=10) :: catime
 namelist/setup/catime,sfcio,lputsi,dsearch,nx,ny

 call init_setup
!
!Read namelist
!
 read(5,setup)

 read(catime(1:10),'(i4,3i2)') iy,im,id,ih

 write(*,*) 'catime,nx,ny = ',catime,nx,ny

!
! force ih = 0 to make it has the same time as th eold SST files
!
 ih = 0

 allocate( tf(nx,ny),sice(nx,ny),tf_clm(nx,ny),sal_clm(nx,ny),mask(nx,ny),masks(nx,ny) )
 allocate( mask_ij(nx*ny),tf_ij(nx*ny))
 allocate( xlats(ny),xlons(nx),xlats_ij(nx*ny),xlons_ij(nx*ny) )
!
! get xlats & xlons
!
  dres = 180.0/real(ny)
  y0 = 0.5*dres-90.0
  x0 = 0.5*dres

! Get lat_sst & lon_sst
  do j = 1, ny
     xlats(j) = y0 + real(j-1)*dres
  enddo

  do i = 1, nx
     xlons(i) = (x0 + real(i-1)*dres)
  enddo
!
! Get 1-dimentional xlats and xlons for global grids (nx*ny)
!
 ij = 0
 do j = 1, ny
    do i = 1, nx
       ij = ij + 1
       xlats_ij(ij) = xlats(j)
       xlons_ij(ij) = xlons(i)
    enddo
 enddo
!
! read target mask (nx,ny): mask = 0 for ocean; mask = 1 for land
!
 call get_mask(fin_grbmsk,mask,nx,ny)
 write(*,*) 'get mask from ',fin_grbmsk

!
! read 0.5 degree mask (nx_0p5,ny_0p5) while processing 1/12 degree
!
 if ( nx == 4320 ) then
 call get_mask(fin_grbmsk_0p5,mask_0p5,nx_0p5,ny_0p5)
    write(*,*) 'get 0.5 degree mask from ',fin_grbmsk_0p5
 endif
!
! Get dimensions of sfcanl, allocate arrays and read NSST tf
!
 if ( trim(sfcio) == 'nemsio' ) then
    call get_sfcanl_dim_nems(fin_sfcanl,ny_nst,nx_nst)
 elseif ( trim(sfcio) == 'ncio' ) then
    call get_sfcanl_dim_nc(fin_sfcanl,ny_nst,nx_nst)
 else
    write(*,*) 'invalid sfcio. abort'
    stop
 endif

 write(*,'(a,8I5)') 'iy,im,id,ih,ny_nst,nx_nst,ny,nx : ',iy,im,id,ih,ny_nst,nx_nst,ny,nx
 allocate( tf_nst(nx_nst,ny_nst),mask_nst(nx_nst,ny_nst) )
 allocate( xlons_nst(nx_nst),xlats_nst(ny_nst) )
 if ( trim(sfcio) == 'nemsio' ) then
    call read_tf_sfcanl_nems(fin_sfcanl,tf_nst,mask_nst,xlats_nst,xlons_nst,ny_nst,nx_nst)
 elseif ( trim(sfcio) == 'ncio' ) then
    call read_tf_sfcanl_nc(fin_sfcanl,tf_nst,mask_nst,xlats_nst,xlons_nst,ny_nst,nx_nst)
 endif
!
! get sea ice analysis (daily and available 6-houly files) at taget grids (nx,ny)
!
  call get_seaice(xlats_ij,xlons_ij,ny,nx,sice)
!
! modify mask for sea ice
!
 masks = mask
 do j = 1, ny
    do i = 1, nx
       if ( sice(i,j) <= 1.0 .and. sice(i,j) >= 0.15 ) then
          mask(i,j) = 2
       endif
    enddo
 enddo
 lensfc=nx*ny  
 mask_ij(:)  = reshape( mask,  (/lensfc/) )
!
! get tf_clm at the analysis time & at target resolution from monthly climatology (nx,ny)
!
  call get_tf_clm(mask_ij,xlats_ij,xlons_ij,ny,nx,iy,im,id,ih,dsearch,miss_fill,bmiss,tf_clm)
!
! get sal_clm at the analysis time & at target resolution from monthly climatology (nx,ny)
!
 call get_sal_clm(xlats_ij,xlons_ij,ny,nx,iy,im,id,ih,sal_clm)
!
! interpolate tf_nst (Gaussian lat/lon) to tf (lat/lon)
!
 call lalo_to_tile(tf_nst, mask_nst,xlats_nst,xlons_nst,ny_nst,nx_nst, &
                   tf_ij,  mask_ij, xlats_ij, xlons_ij, ny,nx, &
                   sfcflag,dsearch, miss_fill,bmiss)
 write(*,*) 'done with lalo_to_tile for tf_ij'

 tf(:,:) = reshape (tf_ij, (/nx,ny/) )
!
! set temp to be RTG SST climatology for the missing water grids 
!
 nfill_clm_water=0
 nfill_clm_land=0
 nfill_clm_ice=0
 do j = 1, ny
    do i = 1, nx
       if ( (mask(i,j) == 0 .and. tf(i,j) == bmiss) .or. mask(i,j) /= 0 ) then
          tf(i,j) = tf_clm(i,j)
!
!         get the counts statistics for a surface type
!
          if ( mask(i,j) == 0 ) then
             nfill_clm_water = nfill_clm_water + 1
          elseif ( mask(i,j) == 1 ) then
             nfill_clm_land = nfill_clm_land + 1
          else
             nfill_clm_ice = nfill_clm_ice + 1
          endif
       endif
    enddo
 enddo

 write(*,*) 'nfill_clm_water : ',nfill_clm_water
 write(*,*) 'nfill_clm_land : ',nfill_clm_land
 write(*,*) 'nfill_clm_ice : ',nfill_clm_ice
!
! calculate water temperature for grids with positive sea ice fraction
! set temperature of some fresh water lakes to 0C
!
 do j = 1, ny
    do i = 1, nx
       if ( mask(i,j) == 2 ) then
          tf(i,j) = tfreez(sal_clm(i,j))+tzero
       endif
    enddo
 enddo

!
! fill land grids temperature with known water grids temperatures
!
 if ( lputsi ) then
    tf(:,:) = tf(:,:) - tzero
    call putsi(tf,mask,nx,ny)
    tf(:,:) = tf(:,:) + tzero
 endif
!
! write out tf as grib1 and grib2 format
!
 maskflag=0      ! no bitmap
 call togrib(fout_tf_grb,tf,masks,maskflag,nx,ny,iy,im,id,ih)

 if ( nx == 4320 ) then
!
! Get half degree Tf
!
    do ii = 6, nx, 6
       i = ii/6
       do jj = 6, ny, 6
          j = jj/6
          tf_0p5(i,j) = 0.0
          do ix = 0, 5
             do jx = 0, 5
                tf_0p5(i,j) = tf_0p5(i,j) + tf(ii-ix,jj-jx)
             enddo
          enddo
          tf_0p5(i,j) = tf_0p5(i,j) / 36.0
       enddo
    enddo
    write(*,*) 'generate half degree grib SST file while processing 1/12 degree file'
    call togrib(fout_tf_grb_0p5,tf_0p5,mask_0p5,maskflag,nx_0p5,ny_0p5,iy,im,id,ih)
!
! Get one degree Tf (for gdas_oisst)
!
    do ii = 12, nx, 12
       i = ii/12
       do jj = 12, ny, 12
          j = jj/12
          tf_1p0(i,j) = 0.0
          rmask_1p0(i,j) = 0.0
          do ix = 0, 11
             do jx = 0, 11
                tf_1p0(i,j) = tf_1p0(i,j) + tf(ii-ix,jj-jx)
                rmask_1p0(i,j) = rmask_1p0(i,j) + mask(ii-ix,jj-jx)
             enddo
          enddo
          tf_1p0(i,j) = tf_1p0(i,j) / 144.0
          rmask_1p0(i,j) = rmask_1p0(i,j) / 144.0
          mask_1p0(i,j) = nint(rmask_1p0(i,j))
       enddo
    enddo
    write(*,*) 'generate one degree grib SST file while processing 1/12 degree file'
    call togrib(fout_tf_grb_1p0,tf_1p0,mask_1p0,maskflag,nx_1p0,ny_1p0,iy,im,id,ih)
 endif

 maskflag=1      ! with bitmap
 if ( nx /= 1440 ) then
    call togrib(fout_tf_grb_awips,tf,masks,maskflag,nx,ny,iy,im,id,ih)
    if ( nx == 4320 ) then
       write(*,*) 'generate half degree bitmap grib SST file while processing 1/12 degree file'
       call togrib(fout_tf_grb_0p5_awips,tf_0p5,mask_0p5,maskflag,nx_0p5,ny_0p5,iy,im,id,ih)
    endif
 endif

 write(*,*) 'All Done'

end program nsst 

