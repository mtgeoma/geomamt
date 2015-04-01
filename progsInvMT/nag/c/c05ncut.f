      SUBROUTINE C05NCU(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14C REVISED. IER-869 (NOV 1990).
C     **********
C
C     SUBROUTINE C05NCU (based on MINPACK routine DOGLEG)
C
C     Given an M by N matrix A, an N by N nonsingular diagonal
C     matrix D, an M-vector B, and A positive number DELTA, the
C     problem is to determine the convex combination X of the
C     Gauss-Newton and scaled gradient directions that minimizes
C     (A*X - B) in the least squares sense, subject to the
C     restriction that the euclidean norm of D*X be at most DELTA.
C
C     This subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     QR factorization of A. That is, if A = Q*R, where Q has
C     orthogonal columns and R is an upper triangular matrix,
C     then C05NCU expects the full upper triangle of R and
C     the first N components of (Q transpose)*B.
C
C     The subroutine statement is
C
C     SUBROUTINE C05NCU(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C
C     where
C
C     N is a positive integer input variable set to the order of R.
C
C     R is an input array of length LR which must contain the upper
C     triangular matrix R stored by rows.
C
C     LR is a positive integer input variable not less than
C     (N*(N+1))/2.
C
C     DIAG is an input array of length N which must contain the
C     diagonal elements of the matrix D.
C
C     QTB is an input array of length N which must contain the first
C     N elements of the vector (Q transpose)*B.
C
C     DELTA is a positive input variable which specifies an upper
C     bound on the Euclidean norm of D*X.
C
C     X is an output array of length N which contains the desired
C     convex combination of the Gauss-Newton direction and the
C     scaled gradient direction.
C
C     WA1 and WA2 are work arrays of length N.
C
C     Argonne National Laboratory. MINPACK project. March 1980.
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
      DOUBLE PRECISION  DELTA
      INTEGER           LR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DIAG(N), QTB(N), R(LR), WA1(N), WA2(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BNORM, EPSMCH, GNORM, QNORM, SGNORM, TEMP
      INTEGER           I, J, JJ, L
C     .. External Functions ..
      DOUBLE PRECISION  F06EJF, X02AJF
      EXTERNAL          F06EJF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          DTPMV, DTPSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
      EPSMCH = X02AJF()
C
C     First, calculate the Gauss-Newton direction.
C
      JJ = 1
      DO 40 J = 1, N
         WA1(J) = R(JJ)
         IF (R(JJ).EQ.ZERO) THEN
            TEMP = ZERO
            L = J
            DO 20 I = 1, J - 1
               TEMP = MAX(TEMP,ABS(R(L)))
               L = L + N - I
   20       CONTINUE
            IF (TEMP.EQ.ZERO) THEN
               R(JJ) = EPSMCH
            ELSE
               R(JJ) = EPSMCH*TEMP
            END IF
         END IF
         JJ = JJ + N - J + 1
   40 CONTINUE
      DO 60 I = 1, N
         X(I) = QTB(I)
   60 CONTINUE
      CALL DTPSV('Lower triangle','Transpose','Non-unit diagonal',N,R,X,
     *           1)
      JJ = 1
      DO 80 J = 1, N
         R(JJ) = WA1(J)
         JJ = JJ + N - J + 1
   80 CONTINUE
C
C     Test whether the Gauss-Newton direction is acceptable.
C
      DO 100 J = 1, N
         WA2(J) = DIAG(J)*X(J)
  100 CONTINUE
      QNORM = F06EJF(N,WA2,1)
      IF (QNORM.GT.DELTA) THEN
C
C        The Gauss-Newton direction is not acceptable.
C        next, calculate the scaled gradient direction.
C
         DO 120 I = 1, N
            WA1(I) = QTB(I)
  120    CONTINUE
         CALL DTPMV('Lower triangle','Not transpose',
     *              'Non-unit diagonal',N,R,WA1,1)
         DO 140 I = 1, N
            WA1(I) = WA1(I)/DIAG(I)
  140    CONTINUE
C
C        Calculate the norm of the scaled gradient and test for
C        the special case in which the scaled gradient is zero.
C
         GNORM = F06EJF(N,WA1,1)
         SGNORM = ZERO
         ALPHA = DELTA/QNORM
         IF (GNORM.NE.ZERO) THEN
C
C           Calculate the point along the scaled gradient
C           at which the quadratic is minimized.
C
            DO 160 J = 1, N
               WA1(J) = (WA1(J)/GNORM)/DIAG(J)
  160       CONTINUE
            DO 180 I = 1, N
               WA2(I) = WA1(I)
  180       CONTINUE
            CALL DTPMV('Lower triangle','Transpose','Non-unit diagonal',
     *                 N,R,WA2,1)
            TEMP = F06EJF(N,WA2,1)
            SGNORM = (GNORM/TEMP)/TEMP
C
C           Test whether the scaled gradient direction is acceptable.
C
            ALPHA = ZERO
            IF (SGNORM.LT.DELTA) THEN
C
C              The scaled gradient direction is not acceptable.
C              Finally, calculate the point along the dogleg
C              at which the quadratic is minimized.
C
               BNORM = F06EJF(N,QTB,1)
               TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
               TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2 +
     *                SQRT((TEMP-(DELTA/QNORM))**2+(ONE-(DELTA/QNORM)
     *                **2)*(ONE-(SGNORM/DELTA)**2))
               ALPHA = ((DELTA/QNORM)*(ONE-(SGNORM/DELTA)**2))/TEMP
            END IF
         END IF
C
C        Form appropriate convex combination of the Gauss-Newton
C        direction and the scaled gradient direction.
C
         TEMP = (ONE-ALPHA)*MIN(SGNORM,DELTA)
         DO 200 J = 1, N
            X(J) = TEMP*WA1(J) + ALPHA*X(J)
  200    CONTINUE
      END IF
      RETURN
      END
