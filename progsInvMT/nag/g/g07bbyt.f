      SUBROUTINE G07BBY(XMU,XSIG,N,TOLM,MAXITS,SUM,SUM2,L11,L12,L22,K,
     *                  WK,IERROR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     ROUTINE FOR COMPUTING MAXIMUM-LIKELIHOOD ESTIMATES
C     FOR PARAMETERS OF THE NORMAL DISTRIBUTION FROM
C     GROUPED AND CENSORED DATA USING THE EM ALGORITHM.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  L11, L12, L22, SUM, SUM2, TOLM, XMU, XSIG
      INTEGER           IERROR, K, MAXITS
C     .. Array Arguments ..
      DOUBLE PRECISION  WK(*)
      INTEGER           N(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, F, L1, L2, RN, RP, SUMG2, TD, TS, W, XNEW,
     *                  YD, YN, YNEW, YNU, YS, YSU
      INTEGER           I, IFAIL, N1, N2, N3
C     .. External Functions ..
      DOUBLE PRECISION  G01MAZ, G01MBF, S15ABF
      EXTERNAL          G01MAZ, G01MBF, S15ABF
C     .. External Subroutines ..
      EXTERNAL          G07BBX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Executable Statements ..
C
      N1 = N(1)
      N2 = N1 + N(2)
      N3 = N2 + N(3) + N(3)
      IFAIL = 0
      IERROR = 0
      K = 0
      RP = DBLE(N(4))
      RN = DBLE(N2+N(3)+N(4))
   20 CONTINUE
C
      TS = SUM
      SUMG2 = SUM2
      TD = RP
      DO 40 I = 1, N1
         YS = (WK(I)-XMU)/XSIG
         F = G01MBF(YS)
         W = XMU + XSIG*F
         TD = TD + F*(F-YS)
         TS = TS + W
         SUMG2 = SUMG2 + (W**2)
   40 CONTINUE
      DO 60 I = N1 + 1, N2
         YS = (WK(I)-XMU)/XSIG
         F = G01MBF(-YS)
         W = XMU - XSIG*F
         TD = TD + F*(F+YS)
         TS = TS + W
         SUMG2 = SUMG2 + (W**2)
   60 CONTINUE
      DO 80 I = N2 + 1, N3, 2
         YS = (WK(I)-XMU)/XSIG
         YSU = (WK(I+1)-XMU)/XSIG
         YN = G01MAZ(YS)
         YNU = G01MAZ(YSU)
         YD = S15ABF(YSU,IFAIL) - S15ABF(YS,IFAIL)
         IF (YD.LE.0.0D0) THEN
            A = 0.0D0
         ELSE
            A = (YN-YNU)/YD
         END IF
         W = XMU + XSIG*A
         TD = TD + ((A**2)+(YSU*YNU-YS*YN)/YD)
         TS = TS + W
         SUMG2 = SUMG2 + (W**2)
   80 CONTINUE
C
      XNEW = TS/RN
      IF (TD.EQ.0.0D0) THEN
         GO TO 120
      ELSE IF ((SUMG2+RN*(XMU**2)-2.0D0*TS*XMU)/TD.LE.0.0D0) THEN
         GO TO 120
      ELSE
         YNEW = SQRT((SUMG2+RN*(XMU**2)-2.0D0*TS*XMU)/TD)
         K = K + 1
         IF (XNEW.EQ.0.0D0) THEN
            IF (ABS(XMU).LT.TOLM .AND. ABS(YNEW-XSIG)/YNEW.LT.TOLM)
     *          GO TO 100
         ELSE IF (ABS((XNEW-XMU)/XNEW).LT.TOLM .AND. ABS((YNEW-XSIG)
     *            /YNEW).LT.TOLM) THEN
            GO TO 100
         END IF
         IF (K.LT.MAXITS) THEN
            XMU = XNEW
            XSIG = YNEW
            GO TO 20
         END IF
      END IF
      IERROR = 2
C
  100 XMU = XNEW
      XSIG = YNEW
      CALL G07BBX(XMU,XSIG,N,SUM,SUM2,L1,L2,L11,L12,L22,WK)
      RETURN
C
  120 IERROR = 3
      RETURN
      END
