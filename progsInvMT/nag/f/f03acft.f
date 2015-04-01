      SUBROUTINE F03ACF(A,IA,N,M,DET,RL,IL,M1,LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 17 REVISED. IER-1673 (JUN 1995).
C
C     DETERMINANT OF REAL POSITIVE DEFINITE SYMMETRIC BAND MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03ACF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DET
      INTEGER           IA, IL, LFAIL, M, M1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,M1), RL(IL,M1)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, XXXX
      INTEGER           ID2, ISAVE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          P01ABF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          F03AGZ
C     .. Intrinsic Functions ..
      INTRINSIC         LOG
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 1
      CALL F03AGZ(N,M,A,IA,RL,IL,M1,D1,ID2,LFAIL)
      IF (LFAIL.EQ.0) GO TO 20
      LFAIL = P01ABF(ISAVE,LFAIL,SRNAME,0,P01REC)
      RETURN
   20 IF (ID2.GT.-LOG(X02AMF())/LOG(2.0D0)) GO TO 40
      IF (ID2.LT.LOG(X02AMF())/LOG(2.0D0)+4) GO TO 60
      DET = D1*2.0D0**ID2
      RETURN
   40 LFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
      RETURN
   60 LFAIL = P01ABF(ISAVE,3,SRNAME,0,P01REC)
      RETURN
      END