      SUBROUTINE G03ACZ(E,LDE,WSUM,NCV,NX,NY,IERROR)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     COMPUTES TEST STATISTICS FOR CANONICAL CORRELATIONS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  WSUM
      INTEGER           IERROR, LDE, NCV, NX, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  E(LDE,6)
C     .. Local Scalars ..
      DOUBLE PRECISION  ADJ, AL, SUML
      INTEGER           I, IFAULT
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      EXTERNAL          G01ECF
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, DBLE
C     .. Executable Statements ..
      IERROR = 0
      SUML = 0.0D0
      DO 20 I = 1, NCV
         AL = E(I,1)*E(I,1)
         IF (AL.GE.1.0D0) THEN
            IERROR = 1
            RETURN
         END IF
         E(I,4) = 1.0D0 - AL
         E(I,2) = AL/(1.0D0-AL)
         SUML = E(I,2) + SUML
   20 CONTINUE
      DO 40 I = 1, NCV
         E(I,3) = E(I,2)/SUML
   40 CONTINUE
      E(NCV,4) = LOG(E(NCV,4))
      DO 60 I = NCV - 1, 1, -1
         E(I,4) = LOG(E(I,4)) + E(I+1,4)
   60 CONTINUE
      ADJ = -(WSUM-DBLE(NX+NY+3)*0.5D0)
      DO 80 I = 1, NCV
         E(I,4) = E(I,4)*ADJ
         E(I,5) = DBLE((NX-I+1)*(NY-I+1))
         IFAULT = 1
         E(I,6) = G01ECF('UPPER',E(I,4),E(I,5),IFAULT)
   80 CONTINUE
      END
