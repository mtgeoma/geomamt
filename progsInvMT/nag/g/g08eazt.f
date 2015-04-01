      SUBROUTINE G08EAZ(NTOT,NR,EX,C,LDC,WRK,LWRK,IERROR)
C     .. Parameters ..
      DOUBLE PRECISION  TWO, ONE
      PARAMETER         (TWO=2.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           IERROR, LDC, LWRK, NR, NTOT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,NR), EX(NR), WRK(LWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  FI1, FJ1, P, RN, S, Y
      INTEGER           I, IF2, J, K
C     .. External Functions ..
      DOUBLE PRECISION  S14AAF
      EXTERNAL          S14AAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      IF (NTOT.LT.1) THEN
         IERROR = 1
      ELSE IF (NR.LT.1) THEN
         IERROR = 1
      ELSE IF (LDC.LT.NR) THEN
         IERROR = 1
      ELSE IF (LWRK.LT.(NR*(NR+5)/2+1)) THEN
         IERROR = 1
      ELSE
C
C        Prepare all the needed factorials
C
         DO 20 I = 1, 2*NR + 1
            IF2 = 1
            WRK(I) = S14AAF(DBLE(I+1),IF2)
   20    CONTINUE
C
C        Calculate the means and the covariance matrix in two stages.
C
C        Stage 1.
C
         RN = DBLE(NTOT+1)
         K = 2*NR + 1
         DO 60 I = 1, NR
            FI1 = WRK(I+1)
            Y = DBLE(I)
            EX(I) = ((RN*Y)+1.0D0-(Y*Y))/FI1
            DO 40 J = 1, I
               K = K + 1
               FJ1 = WRK(J+1)
               IF ((I+J).GT.NTOT) THEN
                  WRK(K) = EX(I) - EX(I)*EX(J)
               ELSE
                  S = DBLE(I+J)
                  P = DBLE(I*J)
                  WRK(K) = EX(I) + RN*(((S*(1-P)+P)/(FI1*FJ1))
     *                     -(TWO*S/WRK(I+J+1))) + (TWO*(S-1)/WRK(I+J)) +
     *                     ((S*S-S-TWO)*P-S*S-P*P+ONE)/(FI1*FJ1)
               END IF
   40       CONTINUE
   60    CONTINUE
C
C        Stage 2.
C
         K = 2*NR + 1
         DO 100 I = 1, NR
            IF (I.LT.NR) EX(I) = EX(I) - EX(I+1)
            DO 80 J = 1, I
               K = K + 1
               IF (I.LT.NR) THEN
                  IF (J.LT.I) THEN
                     C(I,J) = WRK(K) - WRK(K+I) - WRK(K+1) + WRK(K+I+1)
                  ELSE
                     C(I,J) = WRK(K) + WRK(K+I+1) - 2*WRK(K+I)
                  END IF
                  C(J,I) = C(I,J)
               ELSE IF (J.LT.I) THEN
                  C(I,J) = WRK(K) - WRK(K+1)
                  C(J,I) = C(I,J)
               ELSE
                  C(I,J) = WRK(K)
                  C(J,I) = C(I,J)
               END IF
   80       CONTINUE
  100    CONTINUE
      END IF
C
      END
