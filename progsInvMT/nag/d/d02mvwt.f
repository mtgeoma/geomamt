      SUBROUTINE D02MVW(X,XOUT,YOUT,NEQ,KOLD,PHI,NYH,ND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C-----------------------------------------------------------------------
C     THE METHODS IN SUBROUTINE SPDASL USE POLYNOMIALS
C     TO APPROXIMATE THE SOLUTION. D02MVW APPROXIMATES THE
C     solution and its derivative at time xout by evaluating
C     one of these polynomials,and its derivative,there.
C     information defining this polynomial is passed from
C     dastep, so ddatrp cannot be used alone.
C
C     the parameters are%
C     x     the current time in the integration.
C     xout  the time at which the solution is desired
C     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT OR
C           THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT
C           (this is output)
C     neq   number of equations
C     kold  order used on last successful step
C     phi   array of scaled divided differences of y
C           DASSL PHI STARTS AT PHI(1,3) PLEASE NOTE
C     ND    SOLUTION DERIVE REQUIRED MUST BE EITHER 0 OR 1
C
C
C   VIA COMMON....
C     psi   array of past stepsize history
C-----------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, XOUT
      INTEGER           KOLD, ND, NEQ, NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  PHI(NYH,*), YOUT(*)
C     .. Arrays in Common ..
      DOUBLE PRECISION  DUMMY(28), PSI(7)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D, GAMMA, TEMP1
      INTEGER           I, J, KOLDP1
C     .. Common blocks ..
      COMMON            /BD02MV/PSI, DUMMY
      SAVE              /BD02MV/
C     .. Executable Statements ..
      KOLDP1 = KOLD + 1
      TEMP1 = XOUT - X
      DO 20 I = 1, NEQ
         YOUT(I) = PHI(I,1)*(1-ND)
   20 CONTINUE
      C = 1.0D0
      D = 0.0D0
      GAMMA = TEMP1/PSI(1)
      DO 60 J = 2, KOLDP1
         D = D*GAMMA + C/PSI(J-1)
         C = C*GAMMA
         GAMMA = (TEMP1+PSI(J-1))/PSI(J)
         DO 40 I = 1, NEQ
            YOUT(I) = YOUT(I) + (C*(1-ND)+D*ND)*PHI(I,J)
   40    CONTINUE
   60 CONTINUE
      RETURN
C------END OF SUBROUTINE D02MVW------
      END
