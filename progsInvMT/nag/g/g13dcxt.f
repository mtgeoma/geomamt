      SUBROUTINE G13DCX(P,K,PHI,A,IA,R,RR,RI,INTGR,KR,STAT)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     THIS ROUTINE TESTS WHETHER A GIVEN OPERATOR (AR OR MA)
C     IS STATIONARY AND SETS STAT TO .FALSE. IF THIS IS NOT THE
C     CASE
C
C     .. Scalar Arguments ..
      INTEGER           IA, K, KR, P, R
      LOGICAL           STAT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,KR), PHI(K,P*K+1), RI(KR), RR(KR)
      INTEGER           INTGR(KR)
C     .. Local Scalars ..
      INTEGER           I, IFAIL, J, L
C     .. External Subroutines ..
      EXTERNAL          F02AFF
C     .. Executable Statements ..
C
C     INITIALISE A TO THE ZERO MATRIX
C
      DO 40 I = 1, KR
         DO 20 J = 1, KR
            A(I,J) = 0.0D0
   20    CONTINUE
   40 CONTINUE
C
C     INITIALISE NON-ZERO ELEMENTS OF A
C
      DO 100 L = 1, P
         DO 80 I = 1, K
            DO 60 J = 1, K
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
C     CALCULATE EIGENVALUES OF A
C
      IFAIL = 1
      CALL F02AFF(A,IA,KR,RR,RI,INTGR,IFAIL)
      IF (IFAIL.NE.0) THEN
         STAT = .TRUE.
         RETURN
      END IF
C
C     TEST FOR EIGENVALUES ON OR OUTSIDE UNIT CIRCLE
C
      STAT = .TRUE.
      DO 160 I = 1, KR
         IF (RR(I)*RR(I)+RI(I)*RI(I).GE.1.0D0) THEN
            STAT = .FALSE.
            RETURN
         END IF
  160 CONTINUE
C
      RETURN
      END
