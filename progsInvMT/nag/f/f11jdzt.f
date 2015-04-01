      SUBROUTINE F11JDZ(N,NNZ,A,ICOL,ISTR,RDIAG,OMEGA,Y,X)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     F11JDZ solves a system of equations
C
C              Mx = y
C
C     involving the preconditioning matrix:
C
C                       1                   -1           T
C              M =  ---------- ( D + w L ) D  ( D + w L )
C                    w (2 - w)
C
C     corresponding to symmetric successive-over-relaxation on a linear
C     system Ax = b, where A is a sparse symmetric matrix stored in SCS
C     format.
C
C     In the above D is the diagonal part of A, L is the strictly lower
C     triangular part of A, and w is a relaxation parameter.
C
C     Arguments
C     =========
C
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C
C     NNZ    (input) INTEGER
C            On entry, the number of non-zero elements in the lower
C            triangular part of A.
C
C     A      (input) DOUBLE PRECISION array, dimension (NNZ)
C            On entry, the non-zero elements in the lower triangular
C            part of the matrix A, ordered by increasing row index
C            and by increasing column index within each row. Multiple
C            entries for the same row and column indices are not
C            allowed.
C
C     ICOL   (input) INTEGER array, dimension (NNZ)
C            On entry, the column indices corresponding to the non-zero
C            elements given in the array A.
C
C     ISTR   (input) INTEGER array, dimension (N+1)
C            On entry, ISTR(i) gives the starting address in the arrays
C            A, IROW and ICOL, of row i of the matrix A.
C
C     RDIAG  (input) DOUBLE PRECISION array, dimension (N)  -1
C            On entry, the elements of the diagonal matrix D .
C
C     OMEGA  (input) DOUBLE PRECISION
C            On entry, the relaxation parameter w.
C
C     Y      (input) DOUBLE PRECISION array, dimension (N)
C            On entry, the vector y.
C
C     X      (output) DOUBLE PRECISION array, dimension (N)
C            On exit, the vector x.
C
C     ==================================================================
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OMEGA
      INTEGER           N, NNZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NNZ), RDIAG(N), X(N), Y(N)
      INTEGER           ICOL(NNZ), ISTR(N+1)
C     .. Local Scalars ..
      DOUBLE PRECISION  S, SUM
      INTEGER           I, ICJ, J
C     .. External Subroutines ..
      EXTERNAL          F06FDF
C     .. Executable Statements ..
C
C     x = w (2 - w) y.
C
      S = OMEGA*(2.D0-OMEGA)
      CALL F06FDF(N,S,Y,1,X,1)
C
C                  -1
C     x = (D + w L)  x
C
      DO 40 I = 1, N
         SUM = 0.D0
         DO 20 J = ISTR(I), ISTR(I+1) - 2
            SUM = SUM + A(J)*X(ICOL(J))
   20    CONTINUE
         X(I) = RDIAG(I)*(X(I)-OMEGA*SUM)
   40 CONTINUE
C
C     x = D x
C
      DO 60 I = 1, N
         J = ISTR(I+1) - 1
         X(I) = A(J)*X(I)
   60 CONTINUE
C
C                 -T
C     x = (D + w L)  x
C
      DO 100 I = N, 1, -1
         X(I) = RDIAG(I)*X(I)
         DO 80 J = ISTR(I), ISTR(I+1) - 2
            ICJ = ICOL(J)
            X(ICJ) = X(ICJ) - OMEGA*A(J)*X(I)
   80    CONTINUE
  100 CONTINUE
C
      RETURN
      END
