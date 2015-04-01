      SUBROUTINE Y90RPF(UPLO,DIAG,TTYPE,N,KD,A,LDA,AB,LDAB,D,COND,
     *                  SCALE,DETMAN,DETEXP,DIST,SEED)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DETMAN, SCALE
      INTEGER           DETEXP, DIST, KD, LDA, LDAB, N, TTYPE
      CHARACTER         DIAG, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AB(LDAB,*), D(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      INTEGER           DTYPE, I, IFAIL, J, KL, KU
C     .. Local Arrays ..
      DOUBLE PRECISION  ORDIAG(2), RDIAG(2)
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F01ZCF, Y90RAF, Y90SMF
C     .. Intrinsic Functions ..
      INTRINSIC         ANINT
C     .. Executable Statements ..
      IF (KD.LT.0) KD = 0
      IF (KD.GE.N) KD = N - 1
      IF (Y90WAF(UPLO,'L')) THEN
         KL = KD
         KU = 0
      ELSE
         KL = 0
         KU = KD
      END IF
      DTYPE = 2
      CALL Y90RAF('G',DTYPE,TTYPE,N,N,KL,KU,AB,LDAB,D,RDIAG,ORDIAG,COND,
     *            SCALE,DETMAN,DETEXP,DIST,SEED)
      IF (TTYPE.NE.1) THEN
         DO 40 J = 1, N
            DO 20 I = 1, KD + 1
               AB(I,J) = ANINT(AB(I,J))
               IF (I.EQ.KU+1 .AND. AB(I,J).EQ.ZERO) AB(I,J) = ONE
   20       CONTINUE
   40    CONTINUE
         DETMAN = ONE
         DETEXP = 0
         DO 60 I = 1, N
            DETMAN = DETMAN*AB(KU+1,I)
            CALL Y90SMF(DETMAN,DETEXP,4)
   60    CONTINUE
      END IF
      IF (Y90WAF(DIAG,'U')) THEN
         DO 80 I = 1, N
            AB(KU+1,I) = ONE
   80    CONTINUE
         DETMAN = ONE
         DETEXP = 0
      END IF
      IFAIL = 0
      CALL F01ZCF('Unpack',N,N,KL,KU,A,LDA,AB,LDAB,IFAIL)
      RETURN
      END
