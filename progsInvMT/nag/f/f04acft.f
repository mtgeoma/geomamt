      SUBROUTINE F04ACF(A,IA,B,IB,N,M,LR,C,IC,RL,IL,M1,LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 17 REVISED. IER-1674 (JUN 1995).
C     APPROXIMATE SOLUTION OF A SET OF REAL SYMMETRIC POSITIVE
C     DEFINITE BAND LINEAR EQUATIONS WITH MULTIPLE RIGHT
C     HAND SIDES.
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04ACF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IC, IL, LFAIL, LR, M, M1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,M1), B(IB,LR), C(IC,LR), RL(IL,M1)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1
      INTEGER           ID2, ISAVE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03AGZ, F04ALZ
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 1
      CALL F03AGZ(N,M,A,IA,RL,IL,M1,D1,ID2,LFAIL)
      IF (LFAIL.EQ.0) GO TO 20
      LFAIL = P01ABF(ISAVE,LFAIL,SRNAME,0,P01REC)
      RETURN
   20 CALL F04ALZ(N,M,LR,RL,IL,M1,B,IB,C,IC)
      RETURN
      END
