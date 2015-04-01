      SUBROUTINE F04MAX(A,INI,INJ,NZ,N,B,Z)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS SUBROUTINE CALCULATES Z = ( L + I + L**T )*B, WHERE L
C     IS THE STRICTLY LOWER TRIANGULAR PART OF A.
C
C     IF A COMES FROM F01MAF ( MA31A ) THEN IT HAS BEEN SCALED TO
C     HAVE UNIT DIAGONALS AND SO THIS ROUTINE COMPUTES Z = A*B.
C
C     INITIALIZE Z.
C
C     .. Scalar Arguments ..
      INTEGER           N, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NZ), B(N), Z(N)
      INTEGER           INI(NZ), INJ(NZ)
C     .. Local Scalars ..
      INTEGER           I, J, K
C     .. External Subroutines ..
      EXTERNAL          DCOPY
C     .. Executable Statements ..
      CALL DCOPY(N,B,1,Z,1)
C
C     LOOP OVER THE NON-ZEROS IN A.
C
      IF (NZ.LE.0) RETURN
      DO 20 K = 1, NZ
         I = INI(K)
         J = INJ(K)
         Z(I) = Z(I) + A(K)*B(J)
         Z(J) = Z(J) + A(K)*B(I)
   20 CONTINUE
      RETURN
C
C     END OF F04MAX (MA31H ).
C
      END
