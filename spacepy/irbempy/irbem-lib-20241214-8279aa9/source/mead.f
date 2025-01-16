!***************************************************************************************************
! Copyright 2004 S. Bourdarie
!
! This file is part of IRBEM-LIB.
!
!    IRBEM-LIB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    IRBEM-LIB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
!
      SUBROUTINE MEAD(xx, yy, zz, kp, abx, aby, abz)

C CALCULATES INCREMENTAL HIGH ALTITUDE GEOMAGNETIC FIELD IN GAUSS
C xx, yy, zz   SM coordinates in Re
C tilt      in degrees
C abx, aby, abz in gauss

      IMPLICIT NONE

C INPUT:
      REAL*8  xx, yy, zz, tilt
      INTEGER*4 KP
C OUTPUT:
      REAL*8  abx, aby, abz
C LOCAL:
      REAL*8  a(7,4), b(3,4), c(7,4), x, y, z, t
      INTEGER*4 i, j, k, l, m, n 

      COMMON /dip_ang/tilt
C COEFFICIENTS OF MODEL 1 (MEAD-FAIRFIELD 1975)

      DATA ((A(I,J),I=1,7),J=1,4)/17.93D0,-5.79D0,2.98D0,-2.57D0,
     1-0.30D0,-1.47D0,1.05D0,21.79D0,-7.03D0,3.02D0,-2.99D0,-0.62D0,
     2-1.22D0,0.95D0,33.16D0,-6.39D0,4.30D0,-3.25D0,-0.44D0,-1.27D0,
     30.45D0,39.48D0,-2.91D0,5.17D0,-3.86D0,-1.04D0,-1.29D0,-1.14D0/,
     4((B(K,L),K=1,3),L=1,4)/-10.11D0,-1.98D0,0.09D0,-11.84D0,-2.57D0,
     5-0.28D0,-16.54D0,-3.08D0,0.22D0,-19.10D0,-3.50D0,0.23D0/,
     6((C(M,N),M=1,7),N=1,4)/-9.41D0,15.07D0,13.16D0,8.36D0,7.95D0,
     74.55D0,0.51D0,-11.96D0,17.87D0,15.88D0,9.77D0,9.43D0,5.57D0,
     81.53D0,-19.88D0,20.23D0,22.72D0,13.23D0,11.46D0,6.33D0,0.67D0,
     9-22.90D0,22.70D0,26.50D0,15.54D0,11.00D0,7.36D0,1.85D0/

C ROTATE AXES ABOUT Z BY 4 DEGREES
      X = 0.99756405D0 * XX - 0.06975647D0 * YY
      Y = 0.06975647D0 * XX + 0.99756405D0 * YY
C SCALE COORDINATES AND TILT ANGLE
      x = x / 10.0D0
      y = y / 10.0D0
      z = zz / 10.0D0
      T = TILT / 10.0D0
C CALCULATE EXTERNAL FIELD COMPONENTS
      ABX = A(1,KP) * Z + A(2,KP) * X * Z + T * (A(3,KP) + A(4,KP) * X
     &      + A(5,KP) * X * X + A(6,KP) * Y * Y + A(7,KP) * Z * Z)
      ABY = B(1,KP) * Y * Z + T * (B(2,KP) * Y + B(3,KP) * X * Y)
      ABZ = C(1,KP) + C(2,KP) * X + C(3,KP) * X * X + C(4,KP) * Y * Y
     &      + C(5,KP) * Z * Z + T * (C(6,KP) * Z + C(7,KP) * X * Z)

      RETURN
      END
