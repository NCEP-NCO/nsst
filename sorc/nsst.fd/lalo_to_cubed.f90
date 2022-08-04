subroutine lalo_to_tile(tf_lalo,mask_lalo,dlats_lalo,dlons_lalo,jdim_lalo,idim_lalo, &
                        tf_tile,mask_tile,xlats_tile,xlons_tile,jdim_tile,idim_tile, &
                        sfcflag,dsearch,miss_fill,bmiss)
!--------------------------------------------------------------------------------
! abstract: Interpolate lon/lat grid to fv3 native grid (tf_lalo => tf_tile)
!--------------------------------------------------------------------------------

! Input
!
! tf_lalo     : (idim_lalo,idim_lalo) tf at lat/lon regular grid
! mask_lalo   : (idim_lalo,idim_lalo) mask of tf  at lat/lon regular grid
! dlats_lalo  : (jdim_lalo) latitudes along y direction of lat/lon regular grid points
! dlons_lalo  : (idim_lalo) longitudes along x direction of lat/lon regular grid points
! jdim_lalo   : number of y dimension of tf_lalo
! idim_lalo   : number of x dimension of tf_lalo
! mask_tile   : (jdim_tile*idim_tile) mask of tf at cubed sphere grid
! xlats_tile  : (jdim_tile*idim_tile) latitudes of all tile grid points
! xlons_tile  : (jdim_tile*idim_tile) longitudes of all tile grid points
! jdim_tile   : number of y dimension of tf_tile
! idim_tile   : number of x dimension of tf_tile
! sfcflag     : surface flag (mask type) of the target
! dsearch     : maximum search radius in KM
! miss_fill   : modes to fill the grids (when failed in search step)
! bmiss       : missing value or default value
!
! Output
!
! tf_tile     : (jdim_tile*idim_tile) tf at cubed sphere grid

 implicit none

! input/output
 real,    dimension(idim_lalo,jdim_lalo), intent(in)  :: tf_lalo
 integer, dimension(idim_lalo,jdim_lalo), intent(in)  :: mask_lalo
 real,    dimension(jdim_lalo),           intent(in)  :: dlats_lalo
 real,    dimension(idim_lalo),           intent(in)  :: dlons_lalo
 integer, dimension(jdim_tile*idim_tile), intent(in)  :: mask_tile
 real,    dimension(jdim_tile*idim_tile), intent(in)  :: xlats_tile
 real,    dimension(jdim_tile*idim_tile), intent(in)  :: xlons_tile
 integer,                                 intent(in)  :: jdim_lalo,idim_lalo,jdim_tile,idim_tile, &
                                                           sfcflag,miss_fill
 real,                                    intent(in)  :: dsearch,bmiss
 real,    dimension(jdim_tile*idim_tile), intent(out) :: tf_tile

! Local
 real, parameter :: deg2rad=3.1415926/180.0
 real,    dimension(jdim_lalo) :: xlats_lalo
 real,    dimension(idim_lalo) :: xlons_lalo
 real    :: tf,wsum,res_km
 integer :: itile,jtile
 integer :: ii,jj,ij,iii,jjj,nintp,mfill,nfill,nsearch,max_search
 integer :: ilalo,jlalo,ilalop1,jlalop1
 integer :: istart,iend,jstart,jend,krad

 integer, allocatable, dimension(:,:)   :: id1,id2,jdc
 real,    allocatable, dimension(:,:,:) :: agrid,s2c

 print*
 print*,'lalo_to_tile: interpolate from lat/lon grids to any one grid with known lat/lon'

 xlats_lalo = dlats_lalo*deg2rad
 xlons_lalo = dlons_lalo*deg2rad

 allocate(agrid(idim_tile,jdim_tile,2))
 agrid(:,:,1) = reshape (xlons_tile, (/idim_tile,jdim_tile/) )
 agrid(:,:,2) = reshape (xlats_tile, (/idim_tile,jdim_tile/) )
 agrid        = agrid*deg2rad

 allocate(id1(idim_tile,jdim_tile))
 allocate(id2(idim_tile,jdim_tile))
 allocate(jdc(idim_tile,jdim_tile))
 allocate(s2c(idim_tile,jdim_tile,4))

!----------------------------------------------------------------------
! compute bilinear weights for each model point from the nearest
! four lalo points. does not account for mask.  that
! happens later.
!----------------------------------------------------------------------

 call remap_coef( 1, idim_tile, 1, jdim_tile, idim_lalo, jdim_lalo, &
                  xlons_lalo, xlats_lalo, id1, id2, jdc, s2c, agrid )

!----------------------------------------------------------------------
!tf_tile will be output.  initialize to bmiss.
!----------------------------------------------------------------------

 tf_tile = bmiss

 nintp = 0
 nsearch = 0
 nfill   = 0
 mfill   = 0

!----------------------------------------------------------------------
! The maximum distance to search is 500 km. how many gaussian
! grid lengths is that?
!----------------------------------------------------------------------
 res_km = 360.0/real(idim_lalo)*111.0
 max_search  = ceiling(dsearch/res_km)

 write(*,'(a,2F7.1,I8)') 'starting ij_loop,dsearch,res_km,max_search ',dsearch,res_km,max_search

 ij_loop : do ij = 1, jdim_tile*idim_tile

!----------------------------------------------------------------------
! skip non-water points. 
!----------------------------------------------------------------------

    if ( mask_tile(ij) /= 0 ) then
       cycle ij_loop  
    endif

!----------------------------------------------------------------------
!   these are points that are open water 
!----------------------------------------------------------------------
    jtile = (ij-1)/idim_tile + 1
    itile = mod(ij,idim_tile)
    if (itile==0) itile = idim_tile

!----------------------------------------------------------------------
!   see if any of the nearest 4 points mask area open water.  
!   if so, apply tf using bilinear interpolation.
!----------------------------------------------------------------------

    ilalo   = id1(itile,jtile)
    ilalop1 = id2(itile,jtile)
    jlalo   = jdc(itile,jtile)
    jlalop1 = jdc(itile,jtile) + 1

    if ( mask_lalo(ilalo,jlalo)     == sfcflag .or. &
         mask_lalo(ilalop1,jlalo)   == sfcflag .or. &
         mask_lalo(ilalop1,jlalop1) == sfcflag .or. &
         mask_lalo(ilalo,jlalop1)   == sfcflag ) then

       tf = 0.0
       wsum  = 0.0

       if (mask_lalo(ilalo,jlalo) == sfcflag) then
          tf = tf + (s2c(itile,jtile,1)*tf_lalo(ilalo,jlalo))
          wsum  = wsum + s2c(itile,jtile,1)
       endif

       if (mask_lalo(ilalop1,jlalo) == sfcflag) then
          tf = tf + (s2c(itile,jtile,2)*tf_lalo(ilalop1,jlalo))
          wsum  = wsum + s2c(itile,jtile,2)
       endif

       if (mask_lalo(ilalop1,jlalop1) == sfcflag) then
          tf = tf + (s2c(itile,jtile,3)*tf_lalo(ilalop1,jlalop1))
          wsum  = wsum + s2c(itile,jtile,3)
       endif

       if (mask_lalo(ilalo,jlalop1) == sfcflag) then
          tf = tf + (s2c(itile,jtile,4)*tf_lalo(ilalo,jlalop1))
          wsum  = wsum + s2c(itile,jtile,4)
       endif

       if ( wsum > 0.0 ) then
          nintp = nintp + 1
          tf_tile(ij) = tf/wsum
       else
!         write(*,*) 'Warning: No water points nearby, do search/fill'
       endif

    else             ! handle the case when no nearby open water grids 

!----------------------------------------------------------------------
!      no nearby gsi/gaussian open water points. perform a spiral search to
!      find nearest non-land point on gsi/gaussian grid in distance of
!      max_search grid number.
!----------------------------------------------------------------------
       do krad = 1, max_search

          istart = ilalo - krad
          iend   = ilalo + krad
          jstart = jlalo - krad
          jend   = jlalo + krad

          do jj = jstart, jend
             do ii = istart, iend
                if ( (jj == jstart) .or. (jj == jend) .or.   &
                     (ii == istart) .or. (ii == iend))  then
                   if ((jj >= 1) .and. (jj <= jdim_lalo)) then
                      jjj = jj
                      if (ii <= 0) then
                         iii = idim_lalo + ii
                      else if (ii >= (idim_lalo+1)) then
                         iii = ii - idim_lalo
                      else
                         iii = ii
                      end if
                      if (mask_lalo(iii,jjj) == sfcflag) then
                         nsearch = nsearch + 1
                         tf_tile(ij) = tf_lalo(iii,jjj)
                         cycle ij_loop
                      endif ! lalo mask is open water
                   endif
                endif
             enddo
           enddo
        enddo ! krad loop

!       print*,'warning !!!!!! search failed at tile point ',itile,jtile
!
!       Fill with nearby 4-point average of SST climatology
!
!       write(*,*) 'Warning: fill needed'
        if ( miss_fill > 0 ) then 
           nfill = nfill + 1
        else
           mfill = mfill + 1
        endif
     endif
 enddo ij_loop

 write(*,*) 'nintp = ',nintp
 write(*,*) 'nsearch = ',nsearch
 write(*,*) 'nfill   = ',nfill
 write(*,*) 'mfill   = ',mfill
 deallocate(id1, id2, jdc, s2c)

end subroutine lalo_to_tile

subroutine intp_tile(tf_lalo,dlats_lalo,dlons_lalo,jdim_lalo,idim_lalo, &
                     tf_tile,xlats_tile,xlons_tile,jdim_tile,idim_tile)
!--------------------------------------------------------------------------------
! abstract: interpolate lon/lat grid to fv3 native grid (tf_lalo => tf_tile) for
! not dependent on mask
!--------------------------------------------------------------------------------

! input
!
! tf_lalo     : (idim_lalo,idim_lalo) tf at lat/lon regular grid
! dlats_lalo  : (jdim_lalo) latitudes along y direction of lat/lon regular grid
! points
! dlons_lalo  : (idim_lalo) longitudes along x direction of lat/lon regular grid
! points
! jdim_lalo   : number of y dimension of tf_lalo
! idim_lalo   : number of x dimension of tf_lalo
! xlats_tile  : (jdim_tile*idim_tile) latitudes of all tile grid points
! xlons_tile  : (jdim_tile*idim_tile) longitudes of all tile grid points
! jdim_tile   : number of y dimension of tf_tile
! idim_tile   : number of x dimension of tf_tile
!
! output
!
! tf_tile     : (jdim_tile*idim_tile) tf at cubed sphere grid

 implicit none

! input/output
 real,    dimension(idim_lalo,jdim_lalo), intent(in)  :: tf_lalo
 real,    dimension(jdim_lalo),           intent(in)  :: dlats_lalo
 real,    dimension(idim_lalo),           intent(in)  :: dlons_lalo
 real,    dimension(jdim_tile*idim_tile), intent(in)  :: xlats_tile
 real,    dimension(jdim_tile*idim_tile), intent(in)  :: xlons_tile
 integer,                                 intent(in)  :: jdim_lalo,idim_lalo,jdim_tile,idim_tile
 real,    dimension(jdim_tile*idim_tile), intent(out) :: tf_tile

! local
 real, parameter :: deg2rad=3.1415926/180.0
 real,    dimension(jdim_lalo) :: xlats_lalo
 real,    dimension(idim_lalo) :: xlons_lalo
 real    :: tf,wsum,res_km
 integer :: itile,jtile
 integer :: ii,jj,ij,iii,jjj
 integer :: ilalo,jlalo,ilalop1,jlalop1
 integer :: istart,iend,jstart,jend,krad

 integer, allocatable, dimension(:,:)   :: id1,id2,jdc
 real,    allocatable, dimension(:,:,:) :: agrid,s2c

 print*
 print*,'intp_tile: interpolate from lat/lon grids to any one grid with known lat/lon'

 xlats_lalo = dlats_lalo*deg2rad
 xlons_lalo = dlons_lalo*deg2rad

 allocate(agrid(idim_tile,jdim_tile,2))
 agrid(:,:,1) = reshape (xlons_tile, (/idim_tile,jdim_tile/) )
 agrid(:,:,2) = reshape (xlats_tile, (/idim_tile,jdim_tile/) )
 agrid        = agrid*deg2rad

 allocate(id1(idim_tile,jdim_tile))
 allocate(id2(idim_tile,jdim_tile))
 allocate(jdc(idim_tile,jdim_tile))
 allocate(s2c(idim_tile,jdim_tile,4))

!----------------------------------------------------------------------
! compute bilinear weights for each model point from the nearest
! four lalo points. does not account for mask.  that
! happens later.
!----------------------------------------------------------------------
 call remap_coef( 1, idim_tile, 1, jdim_tile, idim_lalo, jdim_lalo, &
                  xlons_lalo, xlats_lalo, id1, id2, jdc, s2c, agrid )

 do ij = 1, jdim_tile*idim_tile

    jtile = (ij-1)/idim_tile + 1
    itile = mod(ij,idim_tile)
    if (itile==0) itile = idim_tile

    ilalo   = id1(itile,jtile)
    ilalop1 = id2(itile,jtile)
    jlalo   = jdc(itile,jtile)
    jlalop1 = jdc(itile,jtile) + 1

    wsum = s2c(itile,jtile,1) + s2c(itile,jtile,2) + &
           s2c(itile,jtile,3) + s2c(itile,jtile,4)

    tf_tile(ij)  = ( s2c(itile,jtile,1)*tf_lalo(ilalo,jlalo)     + &
                     s2c(itile,jtile,2)*tf_lalo(ilalop1,jlalo)   + &
                     s2c(itile,jtile,3)*tf_lalo(ilalop1,jlalop1) + &
                     s2c(itile,jtile,4)*tf_lalo(ilalo,jlalop1) )/wsum
 enddo

 deallocate(id1, id2, jdc, s2c)

end subroutine intp_tile

 subroutine remap_coef( is, ie, js, je,&
      im, jm, lon, lat, id1, id2, jdc, s2c, agrid )

!----------------------------------------------------------------------
! this routine was taken from the forecast model -
! ./atmos_cubed_sphere/tools/fv_treat_da_inc.f90.
!----------------------------------------------------------------------

    implicit none
    integer, intent(in)  :: is, ie, js, je
    integer, intent(in)  :: im, jm
    real,    intent(in)  :: lon(im), lat(jm)
    real,    intent(out) :: s2c(is:ie,js:je,4)
    integer, intent(out), dimension(is:ie,js:je):: id1, id2, jdc
    real,    intent(in)  :: agrid(is:ie,js:je,2)
    ! local:
    real :: rdlon(im)
    real :: rdlat(jm)
    real :: a1, b1
    real, parameter :: pi = 3.1415926
    integer i,j, i1, i2, jc, i0, j0
    do i=1,im-1
      rdlon(i) = 1. / (lon(i+1) - lon(i))
    enddo
    rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

    do j=1,jm-1
      rdlat(j) = 1. / (lat(j+1) - lat(j))
    enddo

    ! * Interpolate to cubed sphere cell center
    do 5000 j=js,je

      do i=is,ie

        if ( agrid(i,j,1)>lon(im) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
        elseif ( agrid(i,j,1)<lon(1) ) then
          i1 = im;     i2 = 1
          a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
        else
          do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
              i1 = i0;  i2 = i0+1
              a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
              go to 111
            endif
          enddo
        endif
111     continue

        if ( agrid(i,j,2)<lat(1) ) then
          jc = 1
          b1 = 0.
        elseif ( agrid(i,j,2)>lat(jm) ) then
          jc = jm-1
          b1 = 1.
        else
          do j0=1,jm-1
            if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
              jc = j0
              b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
              go to 222
            endif
          enddo
        endif
222     continue

        if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
             write(*,*) 'gid=', i,j,a1, b1
        endif

        s2c(i,j,1) = (1.-a1) * (1.-b1)
        s2c(i,j,2) =     a1  * (1.-b1)
        s2c(i,j,3) =     a1  *     b1
        s2c(i,j,4) = (1.-a1) *     b1
        id1(i,j) = i1
        id2(i,j) = i2
        jdc(i,j) = jc
      enddo   !i-loop
5000 continue   ! j-loop

 end subroutine remap_coef

