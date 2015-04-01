      SUBROUTINE E04VDN(MODE,N,NACTIV,NCOLZ,NFREE,NQ,UNITQ,KACTIV,KFREE,
     *                  V,ZY,WRK)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C *********************************************************************
C     E04VDN TRANSFORMS THE VECTOR  V  IN VARIOUS WAYS USING THE
C     MATRIX  Q = ( Z  Y )  DEFINED BY THE INPUT PARAMETERS.
C
C     MODE               RESULT
C     ----               ------
C
C       1                V = Z*V
C       2                V = Y*V
C       3                V = Q*V       (NOT YET USED)
C
C     ON INPUT,  V  IS ASSUMED TO BE ORDERED AS  ( V(FREE)  V(FIXED) ).
C     ON OUTPUT, V  IS A FULL N-VECTOR.
C
C
C       4                V = Z(T)*V
C       5                V = Y(T)*V
C       6                V = Q(T)*V
C
C     ON INPUT,  V  IS A FULL N-VECTOR.
C     ON OUTPUT, V  IS ORDERED AS  ( V(FREE)  V(FIXED) ).
C
C       7                V = Y(T)*V
C       8                V = Q(T)*V
C
C     ON INPUT,  V  IS A FULL N-VECTOR.
C     ON OUTPUT, V  IS AS IN MODES 5 AND 6 EXCEPT THAT V(FIXED)
C     IS NOT SET.
C
C     BEWARE THAT  NCOLZ  WILL SOMETIMES BE  NCOLR.
C     ALSO, MODES  1, 4, 7 AND 8  DO NOT INVOLVE  V(FIXED).
C     NACTIV  AND  THE ARRAY  KACTIV  ARE NOT USED FOR THOSE CASES.
C     ORIGINAL VERSION  APRIL 1983. MODES 7 AND 8 ADDED  APRIL 1984.
C *********************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           MODE, N, NACTIV, NCOLZ, NFREE, NQ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  V(N), WRK(N), ZY(NQ,NQ)
      INTEGER           KACTIV(N), KFREE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
      INTEGER           J, J1, J2, K, KA, KW, L, LENV, NFIXED
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FBF, DAXPY, DCOPY
C     .. Data statements ..
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
      J1 = 1
      J2 = NFREE
      IF (MODE.EQ.1 .OR. MODE.EQ.4) J2 = NCOLZ
      IF (MODE.EQ.2 .OR. MODE.EQ.5 .OR. MODE.EQ.7) J1 = NCOLZ + 1
      LENV = J2 - J1 + 1
      IF (MODE.GE.4) GO TO 140
C
C ---------------------------------------------------------------------
C     MODE = 1, 2  OR  3.
C ---------------------------------------------------------------------
      IF (NFREE.GT.0) CALL F06FBF(NFREE,ZERO,WRK,1)
C
C     COPY  V(FIXED)  INTO THE END OF  WRK.
C
      IF (MODE.EQ.1 .OR. NFIXED.EQ.0) GO TO 20
      CALL DCOPY(NFIXED,V(NFREE+1),1,WRK(NFREE+1),1)
C
C     SET  WRK  =  RELEVANT PART OF  ZY * V.
C
   20 IF (LENV.LE.0) GO TO 60
      IF (UNITQ) CALL DCOPY(LENV,V(J1),1,WRK(J1),1)
      IF (UNITQ) GO TO 60
      DO 40 J = J1, J2
         IF (V(J).NE.ZERO) CALL DAXPY(NFREE,V(J),ZY(1,J),1,WRK,1)
   40 CONTINUE
C
C     EXPAND  WRK  INTO  V  AS A FULL N-VECTOR.
C
   60 CALL F06FBF(N,ZERO,V,1)
      IF (NFREE.EQ.0) GO TO 100
      DO 80 K = 1, NFREE
         J = KFREE(K)
         V(J) = WRK(K)
   80 CONTINUE
C
C     COPY  WRK(FIXED)  INTO THE APPROPRIATE PARTS OF  V.
C
  100 IF (MODE.EQ.1 .OR. NFIXED.EQ.0) GO TO 260
      DO 120 L = 1, NFIXED
         KW = NFREE + L
         KA = NACTIV + L
         J = KACTIV(KA)
         V(J) = WRK(KW)
  120 CONTINUE
      GO TO 260
C
C ---------------------------------------------------------------------
C     MODE = 4, 5, 6, 7  OR  8.
C ---------------------------------------------------------------------
C     PUT THE FIXED COMPONENTS OF  V  INTO THE END OF  WRK.
C
  140 IF (MODE.EQ.4 .OR. MODE.GT.6 .OR. NFIXED.EQ.0) GO TO 180
      DO 160 L = 1, NFIXED
         KW = NFREE + L
         KA = NACTIV + L
         J = KACTIV(KA)
         WRK(KW) = V(J)
  160 CONTINUE
C
C     PUT THE FREE  COMPONENTS OF  V  INTO THE BEGINNING OF  WRK.
C
  180 IF (NFREE.EQ.0) GO TO 240
      DO 200 K = 1, NFREE
         J = KFREE(K)
         WRK(K) = V(J)
  200 CONTINUE
C
C     SET  V  =  RELEVANT PART OF  ZY(T) * WRK.
C
      IF (LENV.LE.0) GO TO 240
      IF (UNITQ) CALL DCOPY(LENV,WRK(J1),1,V(J1),1)
      IF (UNITQ) GO TO 240
      DO 220 J = J1, J2
         V(J) = DDOT(NFREE,ZY(1,J),1,WRK,1)
  220 CONTINUE
C
C     COPY THE FIXED COMPONENTS OF  WRK  INTO THE END OF  V.
C
  240 IF (MODE.EQ.4 .OR. MODE.GT.6 .OR. NFIXED.EQ.0) GO TO 260
      CALL DCOPY(NFIXED,WRK(NFREE+1),1,V(NFREE+1),1)
C
  260 RETURN
C
C     END OF E04VDN (ZYPROD)
      END
