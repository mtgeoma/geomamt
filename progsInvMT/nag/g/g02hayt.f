      SUBROUTINE G02HAY(INDW,N,M,X,IX,MM,CUCV,WGT,MAXIT,NITMON,TOL,NIT,
     *                  SA,SZ,SC2,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     WEIGHTS FOR REGRESSION STANDARDIZED WEIGHTS
C
C     .. Parameters ..
      DOUBLE PRECISION  TL, ATL
      PARAMETER         (TL=1.0D-10,ATL=0.9999D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CUCV, TOL
      INTEGER           IFAIL, INDW, IX, M, MAXIT, MM, N, NIT, NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  SA(MM), SC2(MM), SZ(N), WGT(N), X(IX,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  GAM0, SQDEV2, XMD, XME, XSD, ZNR
      INTEGER           I, IFAIL1, IFAIL2, IFUN, J, JJ, L, NFIRST, NN
C     .. External Functions ..
      DOUBLE PRECISION  G02HAV
      EXTERNAL          G02HAV
C     .. External Subroutines ..
      EXTERNAL          G02HAW, G07DAF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
C     PARAMETER CHECK AND INITIALIZATION
C
      IFAIL = 0
      GAM0 = 0.00001D0
      NN = (M+1)*M/2
      IFUN = 2
      IF (INDW.GT.0) IFUN = 1
      NFIRST = N
      IFAIL2 = 0
C
C     INITIAL VALUE FOR SA
C
      DO 20 I = 1, MM
         SA(I) = 0.0D0
   20 CONTINUE
      DO 80 J = 1, M
         CALL G07DAF(NFIRST,X(1,J),WGT,XME,XMD,XSD,IFAIL2)
         SQDEV2 = SQRT(XSD**2+XME**2)
         JJ = (J*J+J)/2
         IF (SQDEV2.GT.TL) GO TO 40
         SA(JJ) = ATL
         GO TO 60
C
   40    SA(JJ) = 1.0D0/SQDEV2
   60    CONTINUE
   80 CONTINUE
C
C     FINAL VALUE FOR SA
C
      IFAIL1 = 0
      CALL G02HAW(X,SA,CUCV,N,M,MM,IX,MAXIT,NITMON,TOL,NIT,SZ,SC2,IFUN,
     *            IFAIL1)
      IF (IFAIL1.EQ.1) THEN
         IFAIL = 1
         RETURN
C
      END IF
C
C     COMPUTE WEIGHTS
C
      DO 100 L = 1, N
         ZNR = SZ(L)
         IF (INDW.GT.0) THEN
            IF (ZNR.LE.GAM0) ZNR = GAM0
            WGT(L) = 1.0D0/ZNR
         END IF
C
         IF (INDW.LT.0) WGT(L) = SQRT(G02HAV(ZNR,IFUN,CUCV))
  100 CONTINUE
      RETURN
C
      END
