      SUBROUTINE G05HDZ(IP,K,PAR,A,RR,RI,INTGR,KR,LIMIT,BIGEST,STAT,
     *                  IFAILA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1686 (JUN 1995).
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGEST, LIMIT
      INTEGER           IFAILA, IP, K, KR
      LOGICAL           STAT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(KR,KR), PAR(IP*K*K), RI(KR), RR(KR)
      INTEGER           INTGR(KR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ROOT
      INTEGER           I, J, KP
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF
      EXTERNAL          A02ABF
C     .. External Subroutines ..
      EXTERNAL          F02AFF, DCOPY, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Initialise A to the zero matrix
C
      DO 20 I = 1, KR
         CALL F06FBF(KR,0.0D0,A(1,I),1)
   20 CONTINUE
C
C     Initialise non-zero elements of A
C
      KP = K*IP
      DO 40 J = 1, K
         CALL DCOPY(KP,PAR(J),K,A(1,J),1)
   40 CONTINUE
C
      DO 80 J = 1, IP - 1
         DO 60 I = 1, K
            A((J-1)*K+I,J*K+I) = 1.0D0
   60    CONTINUE
   80 CONTINUE
C
C     Calculate eigenvalues of A.
C
      IFAILA = 1
      CALL F02AFF(A,KR,KR,RR,RI,INTGR,IFAILA)
      STAT = .TRUE.
      IF (IFAILA.EQ.0) THEN
C
C        Test for eigenvalues on or outside unit circle.
C
         BIGEST = 0.0D0
         DO 100 I = 1, KR
            ROOT = A02ABF(RR(I),RI(I))
            IF (ROOT.GE.LIMIT) THEN
               STAT = .FALSE.
               GO TO 120
            ELSE
               BIGEST = MAX(BIGEST,ROOT)
            END IF
  100    CONTINUE
      END IF
  120 RETURN
C
      END
