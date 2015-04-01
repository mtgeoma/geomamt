      SUBROUTINE G13DCW(A,IK,K,BOOL,L)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     THIS ROUTINE WILL TEST WHETHER A IS POSITIVE DEFINITE
C     BY FACTORISING A AS LL'.
C
C     .. Scalar Arguments ..
      INTEGER           IK, K
      LOGICAL           BOOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IK,K), L(IK,K)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, ID, IFAIL, J
C     .. External Subroutines ..
      EXTERNAL          F03AEF
C     .. Executable Statements ..
      DO 40 I = 1, K
         DO 20 J = I, K
            A(I,J) = A(J,I)
   20    CONTINUE
   40 CONTINUE
C
      BOOL = .TRUE.
      IFAIL = 1
      CALL F03AEF(K,A,IK,L,SUM,ID,IFAIL)
C
      IF (IFAIL.NE.0) THEN
         BOOL = .FALSE.
         DO 80 I = 1, K
            DO 60 J = 1, I
               A(I,J) = A(J,I)
   60       CONTINUE
   80    CONTINUE
         RETURN
      ELSE
         DO 120 I = 1, K
            DO 100 J = I, 1, -1
               IF (I.GT.J) THEN
                  L(I,J) = A(I,J)
                  A(I,J) = A(J,I)
                  L(J,I) = 0.0D0
C
               ELSE
                  L(I,I) = 1.0D0/L(I,1)
               END IF
  100       CONTINUE
  120    CONTINUE
C
      END IF
C
      RETURN
      END
