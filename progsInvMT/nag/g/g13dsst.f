      SUBROUTINE G13DSS(P,K,PHI,A,IA,R,RR,RI,INTGR,KR,STAT,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, K, KR, P, R
      LOGICAL           STAT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,KR), PHI(K,P*K+1), RI(KR), RR(KR)
      INTEGER           INTGR(KR)
C     .. Local Scalars ..
      INTEGER           I, J, L
C     .. External Subroutines ..
      EXTERNAL          F02AFF
C     .. Executable Statements ..
C
C     initialise A to the zero matrix
C
      DO 40 J = 1, KR
         DO 20 I = 1, KR
            A(I,J) = 0.0D0
   20    CONTINUE
   40 CONTINUE
C
C     initialise non-zero elements of A
C
      DO 100 J = 1, K
         DO 80 L = 1, P
            DO 60 I = 1, K
               A((L-1)*K+I,J) = PHI(I,(L-1)*K+J)
   60       CONTINUE
   80    CONTINUE
  100 CONTINUE
C
      IF (R.GT.1) THEN
         DO 140 L = 1, R - 1
            DO 120 I = 1, K
               A((L-1)*K+I,L*K+I) = 1.0D0
  120       CONTINUE
  140    CONTINUE
      END IF
C
C     calculate eigenvalues of A
C
      IFAIL = 1
      CALL F02AFF(A,IA,KR,RR,RI,INTGR,IFAIL)
      STAT = .TRUE.
      IF (IFAIL.NE.0) RETURN
C
C     test for eigenvalues on or outside unit circle
C
      DO 160 I = 1, KR
         IF (RR(I)*RR(I)+RI(I)*RI(I).GE.0.999D0) THEN
            STAT = .FALSE.
            RETURN
         END IF
  160 CONTINUE
      RETURN
C
      END
