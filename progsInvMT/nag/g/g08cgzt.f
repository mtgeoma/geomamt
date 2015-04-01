      DOUBLE PRECISION FUNCTION G08CGZ(CDIST,RINT,PAR,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Provides the probabilities from the appropriate CDF.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 RINT
      INTEGER                          IFAIL
      CHARACTER                        CDIST
C     .. Array Arguments ..
      DOUBLE PRECISION                 PAR(2)
C     .. Local Scalars ..
      DOUBLE PRECISION                 P, PEXP, Q, SR, TOL, Z
      INTEGER                          IFA
C     .. External Functions ..
      DOUBLE PRECISION                 S15ABF, X02AJF, X02AMF
      EXTERNAL                         S15ABF, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL                         S14BAF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG, SQRT
C     .. Executable Statements ..
C
      G08CGZ = 0.0D0
      IF (CDIST.EQ.'N' .OR. CDIST.EQ.'n') THEN
         Z = (RINT-PAR(1))/SQRT(PAR(2))
         IFA = 1
         G08CGZ = S15ABF(Z,IFA)
      ELSE IF (CDIST.EQ.'E' .OR. CDIST.EQ.'e') THEN
         SR = LOG(X02AMF())
         PEXP = PAR(1)*RINT
         IF (-PEXP.GE.SR) THEN
            G08CGZ = 1.0D0 - EXP(-PEXP)
         ELSE
            G08CGZ = 1.0D0
         END IF
      ELSE IF (CDIST.EQ.'C' .OR. CDIST.EQ.'c') THEN
         TOL = X02AJF()*10.0D0
         IFA = 1
         CALL S14BAF(PAR(1)/2.0D0,RINT/2.0D0,TOL,P,Q,IFA)
         G08CGZ = P
      ELSE IF (CDIST.EQ.'G' .OR. CDIST.EQ.'g') THEN
         Z = RINT/PAR(2)
         TOL = X02AJF()*10.0D0
         IFA = 1
         CALL S14BAF(PAR(1),Z,TOL,P,Q,IFA)
         G08CGZ = P
      END IF
      IF (IFA.EQ.3) IFAIL = 1
      RETURN
      END
