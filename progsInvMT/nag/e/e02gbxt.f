      SUBROUTINE E02GBX(KCNT,K,N,ZZ,IZR,MPL1,DD,RR,E,IER,IRR,X,F,RES,
     *                  INDX,PSW,W,IW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8A REVISED. IER-256 (AUG 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-814 (DEC 1989).
C
C
C     ***************
C     CL1 VERSION OF E02GBX.
C     THIS ROUTINE ADMINISTERS THE ADJUSTMENT OF THE
C     Z*D*R   DECOMPOSITION FOR ANY NEW ZERO RESIDUALS.
C     DATA IS PERTURBED AS NECESSARY TO RESOLVE DEGENERACIES.
C     THE DATA CORRESPONDING TO THE ZERO RESIDUALS IS INDEXED
C     IN  INDX(K+1),...,INDX(KCNT).
C
C     W  IS A SCRATCH ARRAY.
C
C     EPS  IS THE SMALLEST POSITIVE NUMBER WHICH
C     SATISFIES  (1.0 + EPS) .GT. 1.0  IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     ***************
C
C     .. Scalar Arguments ..
      INTEGER           IER, IRR, IW, IZR, K, KCNT, MPL1, N
      LOGICAL           PSW
C     .. Array Arguments ..
      DOUBLE PRECISION  DD(N), E(IER,MPL1), F(MPL1), RES(MPL1), RR(IRR),
     *                  W(IW), X(N), ZZ(IZR,N)
      INTEGER           INDX(MPL1)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPS
C     .. Local Scalars ..
      DOUBLE PRECISION  COLNRM, PRJNRM
      INTEGER           I, ISTRT, IX, KP1, TOP
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          E02GBR, E02GBS, E02GBW
C     .. Common blocks ..
      COMMON            /AE02GB/EPS
C     .. Executable Statements ..
      TOP = N + 1
      ISTRT = K + 1
      IF (ISTRT.GT.KCNT) GO TO 60
      DO 40 I = ISTRT, KCNT
         KP1 = K + 1
         IX = INDX(I)
         CALL E02GBS(K,N,ZZ,IZR,N,DD,E(1,IX),W,W(TOP))
         COLNRM = DNRM2(N,E(1,IX),1)
         PRJNRM = DNRM2(N,W,1)
         IF (PRJNRM.LE.EPS*COLNRM) GO TO 20
         INDX(I) = INDX(KP1)
         INDX(KP1) = IX
         CALL E02GBR(K,N,ZZ,IZR,IRR,DD,RR,E(1,IX),W,N)
         K = KP1
         GO TO 40
   20    CONTINUE
         CALL E02GBW(N,E(1,IX),F(IX),X,RES(IX))
   40 CONTINUE
   60 CONTINUE
      RETURN
      END
