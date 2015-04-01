      SUBROUTINE C05NCX(N,A,LDA,RDIAG,ACNORM)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05NCX (Based on MINPACK routine QRFAC )
C
C     This subroutine uses householder transformations to compute a QR
C     factorization of the N by N matrix A. That is, C05NCX determines
C     an orthogonal matrix Q, and an upper triangular matrix R such
C     that A = Q*R. The Householder transformation for column K, K =
C     1,2,...,N, is of the form
C
C                     T
C     I - (1/U(K))*U*U
C
C     where U has zeros in the first K-1 positions. The form of this
C     transformation first appeared in the corresponding LINPACK
C     subroutine.
C
C     The subroutine statement is
C
C        SUBROUTINE C05NCX(N,A,LDA,RDIAG,ACNORM)
C
C     where
C
C     N is a positive integer input variable set to the number
C     of rows and columns of A.
C
C     A is an N by N array. On input A contains the matrix for
C     which the QR factorization is to be computed. On output
C     the strict upper triangular part of a contains the strict
C     upper triangular part of R, and the lower triangular
C     part of a contains a factored form of Q (the non-trivial
C     elements of the U vectors described above).
C
C     LDA is a positive integer input variable not less than N
C     which specifies the leading dimension of the array A.
C
C     RDIAG is an output array of length N which contains the
C     diagonal elements of R.
C
C     ACNORM is an output array of length N which contains the
C     norms of the corresponding columns of the input matrix A.
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
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), ACNORM(N), RDIAG(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJNORM, EPSMCH
      INTEGER           I, J, LDA1
C     .. External Functions ..
      DOUBLE PRECISION  F06EJF, X02AJF
      EXTERNAL          F06EJF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER
C     .. Executable Statements ..
      EPSMCH = X02AJF()
C
C     Compute the intial column norms
C
      DO 20 J = 1, N
         ACNORM(J) = F06EJF(N,A(1,J),1)
   20 CONTINUE
C
C     Reduce A to R with Householder transformations.
C
      LDA1 = LDA
      DO 60 J = 1, N
C
C        Compute the Householder transformation to reduce the
C        J-th column of A to a multiple of the J-th unit vector.
C
         AJNORM = F06EJF(N-J+1,A(J,J),1)
         IF (AJNORM.NE.ZERO) THEN
            IF (A(J,J).LT.ZERO) AJNORM = -AJNORM
            DO 40 I = J, N
               A(I,J) = A(I,J)/AJNORM
   40       CONTINUE
            A(J,J) = A(J,J) + ONE
            IF (J.NE.N) THEN
               IF (J.EQ.N-1) LDA1 = 2
               CALL DGEMV('Transpose',N-J+1,N-J,ONE,A(J,J+1),LDA1,A(J,J)
     *                    ,1,ZERO,RDIAG(J+1),1)
               CALL DGER(N-J+1,N-J,-ONE/A(J,J),A(J,J),1,RDIAG(J+1),1,
     *                   A(J,J+1),LDA1)
            END IF
         END IF
         RDIAG(J) = -AJNORM
   60 CONTINUE
      RETURN
      END
