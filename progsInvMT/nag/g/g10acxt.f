      SUBROUTINE G10ACX(U,TOL,MAXCAL,CRIT,RHO,X,AVH,WWT,N,P,Q,YHAT,C,
     *                  LDC,WK,RES,D,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C      Find RHO for a given CRIT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AVH, CRIT, P, Q, RHO, TOL, U
      INTEGER           IERROR, LDC, MAXCAL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,3), D(17), RES(N), WK(7*(N+2)), WWT(N),
     *                  X(N), YHAT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  FX, XX, YY
      INTEGER           I, IFAIL2, IND1
      LOGICAL           WHILE
C     .. External Functions ..
      DOUBLE PRECISION  G10ACY
      EXTERNAL          G10ACY
C     .. External Subroutines ..
      EXTERNAL          C05AZF
C     .. Executable Statements ..
      XX = 0.0D0
      YY = U
      IFAIL2 = 1
      IND1 = 1
      WHILE = .TRUE.
      I = 0
   20 CONTINUE
      IF (WHILE) THEN
         IF (I.LE.MAXCAL) THEN
            CALL C05AZF(XX,YY,FX,TOL,0,D,IND1,IFAIL2)
            IF (IFAIL2.EQ.1 .AND. IND1.EQ.0) THEN
               IERROR = 4
               GO TO 40
            ELSE IF (IFAIL2.EQ.5) THEN
               IERROR = 5
            END IF
            FX = CRIT - G10ACY(1,XX,X,AVH,WWT,N,P,Q,YHAT,C,LDC,WK,
     *           WK(3*N+7),WK(5*N+11),RES)
            WHILE = IND1 .NE. 0
            I = I + 1
            GO TO 20
         ELSE
            IERROR = 6
         END IF
      END IF
      RHO = (XX+YY)/2.0D0
   40 RETURN
      END
