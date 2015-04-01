      SUBROUTINE C05NCY(M,N,A,LDA,V,W)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05NCY (based on MINPACK routine R1MPYQ)
C
C     Given an M by N matrix A, this subroutine computes A*Q where
C     Q is the product of 2*(N - 1) transformations
C
C        GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
C     eliminate elements in the I-th and N-th planes, respectively.
C     Q itself is not given, rather the information to recover the
C     GV, GW rotations is supplied.
C
C     The subroutine statement is
C
C        SUBROUTINE C05NCY(M,N,A,LDA,V,W)
C
C     where
C
C     M is a positive integer input variable set to the number
C     of rows of A.
C
C     N is a positive integer input variable set to the number
C     of columns of A.
C
C     A is an M by N array. On input A must contain the matrix
C     to be postmultiplied by the orthogonal matrix Q
C     described above. On output A*Q has replaced A.
C
C     LDA is a positive integer input variable not less than M
C     which specifies the leading dimension of the array A.
C
C     V is an input array of length N. V(I) must contain the
C     information necessary to recover the Givens rotation GV(I)
C     described above.
C
C     W is an input array of length N. W(I) must contain the
C     information necessary to recover the Givens rotation GW(I)
C     described above.
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C
C     **********
C
C     Revised to call BLAS.
C     P.J.D. Mayes, J.J. Du Croz, NAG Central Office, September 1987.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), V(N), W(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  COSINE, SINE
      INTEGER           J, NM1, NMJ
C     .. External Subroutines ..
      EXTERNAL          DROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
C     Apply the first set of Givens rotations to A.
C
      NM1 = N - 1
      IF (NM1.GE.1) THEN
         DO 20 NMJ = 1, NM1
            J = N - NMJ
            IF (ABS(V(J)).GT.ONE) THEN
               COSINE = ONE/V(J)
               SINE = SQRT(ONE-COSINE**2)
            ELSE
               SINE = V(J)
               COSINE = SQRT(ONE-SINE**2)
            END IF
            CALL DROT(M,A(1,N),1,A(1,J),1,COSINE,SINE)
   20    CONTINUE
C
C        Apply the second set of Givens rotations to A.
C
         DO 40 J = 1, NM1
            IF (ABS(W(J)).GT.ONE) THEN
               COSINE = ONE/W(J)
               SINE = SQRT(ONE-COSINE**2)
            ELSE
               SINE = W(J)
               COSINE = SQRT(ONE-SINE**2)
            END IF
            CALL DROT(M,A(1,J),1,A(1,N),1,COSINE,SINE)
   40    CONTINUE
      END IF
      RETURN
      END
