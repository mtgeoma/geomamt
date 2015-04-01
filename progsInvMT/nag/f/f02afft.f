      SUBROUTINE F02AFF(A,IA,N,RR,RI,INTGER,LFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     EIGENVALUES OF A REAL UNSYMMETRIC MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AFF')
C     .. Scalar Arguments ..
      INTEGER           IA, LFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), RI(N), RR(N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC
      INTEGER           IB, ISAVE, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF, X02BHF
      EXTERNAL          X02AJF, P01ABF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F01AKF, F01ATF, F02APF
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 1
      ACC = X02AJF()
      IB = X02BHF()
      CALL F01ATF(N,IB,A,IA,K,L,RR)
      CALL F01AKF(N,K,L,A,IA,INTGER)
      CALL F02APF(N,ACC,A,IA,RR,RI,INTGER,LFAIL)
      IF (LFAIL.NE.0) LFAIL = P01ABF(ISAVE,LFAIL,SRNAME,0,P01REC)
      RETURN
      END
