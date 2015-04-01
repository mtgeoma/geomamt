      SUBROUTINE G13DBT(A,IA,B,IB,C,IC,NS,IERR)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        G13DBT MULTIPLIES 2 (NS,NS) MATRICES RESIDENT
C        IN TOP LEFT HAND CORNER OF MATRICES B AND C
C        AND STORES THE RESULT IN TOP LEFT HAND CORNER OF
C        MATRIX A. A,B, AND C MUST BE DISTINCT.
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IC, IERR, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,NS), B(IB,NS), C(IC,NS)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          F01CKZ
C     .. Executable Statements ..
      IERR = 0
      IF (IA.LT.NS) GO TO 60
      IF (IB.LT.NS) GO TO 60
      IF (IC.LT.NS) GO TO 60
      DO 40 J = 1, NS
         DO 20 I = 1, NS
            A(I,J) = 0.0D0
   20    CONTINUE
         CALL F01CKZ(B,IB,NS,NS,C(1,J),1,A(1,J))
   40 CONTINUE
      GO TO 80
   60 IERR = 1
   80 CONTINUE
      RETURN
      END
