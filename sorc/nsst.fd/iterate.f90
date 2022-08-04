      SUBROUTINE iterate(y, mask, nx, ny, dlat, dlon, landflag)
      IMPLICIT none
      INTEGER nx, ny, landflag
      REAL dlat, dlon
      REAL y(nx, ny)
      INTEGER mask(nx, ny)

      INTEGER i, j
      REAL rdlat, rdlon, firstlat
      REAL theta, divisor, del, rpdg, c1, c2, c3
      INTEGER ipi, ipj, imi, imj, jpi, jpj, jmi, jmj

      rpdg = ABS(acos(-1.)/180.)
      rdlat = dlat * rpdg
      rdlon = dlon * rpdg
      c1 = 1./rdlat/rdlat

      firstlat = -90. + dlat/2.

! Need to take care of j = 1, ny seam/edge
!  well, not really -- in antarctic, it's a fool who thinks the pole is
!  not land, and in the arctic, it's analyzed.
      DO j = ny-1, 2, -1
        ipj = j
        jpj = j + 1
        imj = j
        jmj = j - 1
        theta = rpdg * (firstlat + dlat*(j-1))
        divisor = 2./cos(theta)/cos(theta)/rdlon/rdlon + 2./rdlat/rdlat
        c2 = 1./rdlon/rdlon/cos(theta)/cos(theta)
        c3 = tan(theta)/2./rdlat

      DO i = 1, nx

        IF (mask(i,j) .EQ. landflag) THEN
          jmi = i
          jpi = i
          imi = i-1
          IF (imi .EQ. 0) imi = nx
          ipi = i+1
          IF (ipi .EQ. nx+1) ipi = 1
          del = (y(jpi, jpj)+y(jmi, jmj))*c1 + &
                (y(ipi, ipj)+y(imi, imj))*c2
          del = del - c3*(y(jpi, jpj)-y(jmi, jmj))
          y(i,j) = del/divisor
        ENDIF

      ENDDO
      ENDDO


      RETURN
      END
