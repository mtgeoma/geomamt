      SUBROUTINE D05BAQ(WKY,STWT,VK,VG,VF,VW,LSTWT,NSTART,INTIAL,NI,INC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     ------------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutines evaluates the INHOMoGenious terms for
C     n = INTIAL, INTIAL+1, ..., NI.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           INC, INTIAL, LSTWT, NI, NSTART
C     .. Array Arguments ..
      DOUBLE PRECISION  STWT(LSTWT,0:NSTART), VF(0:INC*NI),
     *                  VG(0:NSTART), VK(0:INC*NI), VW(0:NSTART),
     *                  WKY(0:NI)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUMIN
      INTEGER           I, J
C     .. Executable Statements ..
C
      IF (NI.LE.LSTWT) THEN
         DO 40 I = INTIAL, NI
            SUMIN = VF(INC*I)
            DO 20 J = 0, NSTART
               SUMIN = SUMIN + STWT(I,J)*VK(INC*(I-J))*VG(J)
   20       CONTINUE
            WKY(I) = SUMIN
   40    CONTINUE
      ELSE
         IF (INTIAL.LE.LSTWT) THEN
            DO 80 I = INTIAL, LSTWT
               SUMIN = VF(INC*I)
               DO 60 J = 0, NSTART
                  SUMIN = SUMIN + STWT(I,J)*VK(INC*(I-J))*VG(J)
   60          CONTINUE
               WKY(I) = SUMIN
   80       CONTINUE
            DO 100 J = 0, NSTART
               VW(J) = STWT(LSTWT,J)*VG(J)
  100       CONTINUE
            DO 140 I = LSTWT + 1, NI
               SUMIN = VF(INC*I)
               DO 120 J = 0, NSTART
                  SUMIN = SUMIN + VW(J)*VK(INC*(I-J))
  120          CONTINUE
               WKY(I) = SUMIN
  140       CONTINUE
         ELSE
            DO 160 J = 0, NSTART
               VW(J) = STWT(LSTWT,J)*VG(J)
  160       CONTINUE
            DO 200 I = INTIAL, NI
               SUMIN = VF(INC*I)
               DO 180 J = 0, NSTART
                  SUMIN = SUMIN + VW(J)*VK(INC*(I-J))
  180          CONTINUE
               WKY(I) = SUMIN
  200       CONTINUE
         END IF
      END IF
      RETURN
      END
