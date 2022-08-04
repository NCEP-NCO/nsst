      SUBROUTINE sor_lapl(analysis, mask, nx, ny, dlat, dlon, landflag,&
                     weight)
      IMPLICIT none
      INTEGER nx, ny, landflag
      REAL dlat, dlon
      REAL analysis(nx, ny), tmp(nx, ny), old(nx, ny)
      INTEGER mask(nx, ny)
      INTEGER i, itmax
      REAL rms, limit, dmax, dmin, weight

      i      = 0
      itmax  = 20000
      rms    = 9e9
      limit  = 0.001

      old = analysis

      DO WHILE (rms .GT. limit .AND. i .LT. itmax)
        CALL iterate(analysis, mask, nx, ny, dlat, dlon, landflag)
        tmp = analysis - old;

        dmax = MAXVAL(tmp)
        dmin = MINVAL(tmp)
        rms = max(dmax, -dmin);

!        old = analysis;
        analysis = old + weight*tmp; old = analysis
        i = i + 1
!       PRINT *,'starting iteration ',i,'rms = ',rms
      END DO
      PRINT *,'final iteration count, and delta max ',i, rms

      RETURN
      END
