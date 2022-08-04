 real function tfreez(salinity)
!Constants taken from Gill, 1982.
!Author: Robert Grumbine
!LAST MODIFIED: 21 September 1994.

 implicit none

 REAL salinity,sal
 REAL a1, a2, a3
 PARAMETER (a1 = -0.0575)
 PARAMETER (a2 =  1.710523E-3)
 PARAMETER (a3 = -2.154996E-4)

 IF (salinity .LT. 0.) THEN
   sal = 0.
 ELSE
   sal = salinity
 ENDIF
 tfreez = sal*(a1+a2*SQRT(sal)+a3*sal)

 return
 end
