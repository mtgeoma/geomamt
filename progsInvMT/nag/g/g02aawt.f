      SUBROUTINE G02AAW(N,IA,A)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-888 (NOV 1990).
C     .. Scalar Arguments ..
      INTEGER           IA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA*N)
C     .. Local Scalars ..
      INTEGER           I, IERROR, IN, J
C     .. External Subroutines ..
      EXTERNAL          DCOPY
C     .. Executable Statements ..
      IF (N.LE.0) THEN
         IERROR = 1
         RETURN
      ELSE IF (IA.LT.N) THEN
         IERROR = 2
         RETURN
      ELSE IF (N.EQ.2 .AND. IA.EQ.2) THEN
         A(4) = A(3)
         A(3) = A(2)
      ELSE IF (N.GE.2) THEN
         IN = (N*N+N)/2 + 1
         DO 20 I = N, 2, -1
            IN = IN - I
            CALL DCOPY(I,A(IN),1,A(IA*(I-1)+1),1)
   20    CONTINUE
         DO 60 J = 1, N - 1
            DO 40 I = J + 1, N
               A((J-1)*IA+I) = A((I-1)*IA+J)
   40       CONTINUE
   60    CONTINUE
      END IF
      END
