      SUBROUTINE F02AJF(AR,IAR,AI,IAI,N,RR,RI,INTGER,LFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     EIGENVALUES OF A COMPLEX UNSYMMETRIC MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AJF')
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, LFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), RI(N), RR(N)
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
      EXTERNAL          F01AMF, F01AVF, F02ANF
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 1
      ACC = X02AJF()
      IB = X02BHF()
      CALL F01AVF(N,IB,AR,IAR,AI,IAI,K,L,RR)
      CALL F01AMF(N,K,L,AR,IAR,AI,IAI,INTGER)
      CALL F02ANF(N,ACC,AR,IAR,AI,IAI,RR,RI,LFAIL)
      IF (LFAIL.NE.0) LFAIL = P01ABF(ISAVE,LFAIL,SRNAME,0,P01REC)
      RETURN
      END
