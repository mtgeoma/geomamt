      SUBROUTINE G11SAW(Q,XN,AX,LP,IPRINT,ROOTPI)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATE GAUSS-HERMITE QUADRATURE POINTS AND WEIGHTS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ROOTPI
      INTEGER           IPRINT, LP, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(20), XN(20)
C     .. Local Scalars ..
      DOUBLE PRECISION  A2, B2
      INTEGER           IFAIL, ITYPE, K
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          D01BAW, D01BBF, X04BAY
C     .. Executable Statements ..
      A2 = 0.0D0
      B2 = 0.5D0
      ITYPE = 0
      IFAIL = 0
C
      CALL D01BBF(D01BAW,A2,B2,ITYPE,Q,AX,XN,IFAIL)
C
C     CORRECT THE WEIGHTS
C
      DO 20 K = 1, Q
         AX(K) = AX(K)*ROOTPI
   20 CONTINUE
C
      A2 = 0.0D0
      DO 40 K = 1, Q
         A2 = A2 + AX(K)
   40 CONTINUE
      DO 60 K = 1, Q
         AX(K) = AX(K)/A2
   60 CONTINUE
C
      IF (IPRINT.LE.0) RETURN
C
      WRITE (REC,FMT=99999)
      CALL X04BAY(LP,3,REC)
C
      RETURN
C
99999 FORMAT (/' *****************************************************',
     *  '*******************',/)
      END
