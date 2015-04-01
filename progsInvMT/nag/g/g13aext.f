      SUBROUTINE G13AEX(ZB,NZB,EPS,PG,KC)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AEX CHECKS THE VALIDITY OF A SET OF PARAMETER VALUES
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, PG
      INTEGER           KC, NZB
C     .. Array Arguments ..
      DOUBLE PRECISION  ZB(NZB)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ETK, G, Q, RK, U, UME, UPE
      INTEGER           I, J, JH, K, KM, L
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Data statements ..
      DATA              U/1.0D0/
C     .. Executable Statements ..
C
C     CALCULATE TOLERANCE LIMITS ON EITHER SIDE OF UNITY
C
      UPE = U + EPS
      UME = U - EPS
      PG = U
      KC = 1
C
C     A RECURSIVE ALGORITHM IS USED IN WHICH EACH VALUE OF ZB IS
C     EXAMINED IN TURN, WORKING BACKWARDS THROUGH THE SET
C
      DO 140 I = 1, NZB
         K = NZB - I + 1
         KM = K - 1
         A = ZB(K)
C
C        IF THE PARAMETER VALUE EXCEEDS UPE, THE SET IS INVALID
C
         IF (ABS(A).GT.UPE) GO TO 60
C
C        IF THE PARAMETER VALUE IS LESS THAN UME, UPDATE PG AND MODIFY
C        THE REMAINING PARAMETER VALUES
C
         IF (ABS(A).LT.UME) GO TO 80
C
C        THERE IS A ZERO ON THE UNIT CIRCLE, SO WE SET KC TO ZERO
C
         KC = 0
         IF (K.EQ.1) GO TO 160
         RK = DBLE(K)
         ETK = EPS*RK
         DO 20 J = 1, KM
            L = K - J
            Q = ABS(ZB(J)+A*ZB(L))
C
C           IF THE FOLLOWING HOLDS THERE IS AN INVALIDITY CONDITION
C           AND KC IS SET TO -1
C
            IF (Q.GE.ETK) GO TO 60
   20    CONTINUE
C
C        MODIFY THE REMAINING PARAMETER VALUES
C
         DO 40 J = 1, KM
            ZB(J) = DBLE(K-J)*ZB(J)/RK
   40    CONTINUE
         GO TO 140
   60    KC = -1
         GO TO 160
   80    G = (U-A)*(U+A)
         PG = PG*G
         IF (K.EQ.1) GO TO 160
C
C        MODIFY THE REMAINING PARAMETER VALUES
C
         JH = KM/2
         IF (JH.LE.0) GO TO 120
         DO 100 J = 1, JH
            L = K - J
            Q = ZB(J) + A*ZB(L)
            ZB(L) = (ZB(L)+A*ZB(J))/G
            ZB(J) = Q/G
  100    CONTINUE
  120    IF (KM.EQ.(2*JH)) GO TO 140
         ZB(JH+1) = ((U+A)*ZB(JH+1))/G
  140 CONTINUE
  160 RETURN
      END
