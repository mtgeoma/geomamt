      SUBROUTINE F01AEF(N,A,IA,B,IB,DL,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     REDUC1
C     Reduction of the general symmetric eigenvalue problem
C        A * X = LAMBDA * B * X,
C     with symmetric matrix A and symmetric positive definite
C     matrix B, to the equivalent standard problem
C        P * Z = LAMBDA * Z.
C     The upper triangles, including diagonal elements, of A and B
C     are given in the arrays A(N,N) and B(N,N). L(B = L * LT) is
C     formed in the remaining strictly lower triangle of the array
C     B with its diagonal elements in the array DL(N), and the lower
C     triangle of the symmetric matrix P(P = inv(L) * A * inv(LT))
C     is formed in the lower triangle of the array A, including the
C     diagonal elements. Hence the diagonal elements of A are lost.
C     The subroutine will fail if B, perhaps on account of rounding
C     errors, is not positive definite.
C     1st December 1971
C
C     Rewritten to call LAPACK routine SPOTRF/F07FDF;
C     new IFAIL exit inserted for illegal input parameters;
C     error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01AEF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), B(IB,*), DL(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y
      INTEGER           I, IB1, IERR, INFO, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSWAP, DTRSV, F07FDF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
C
      IF (N.LT.0) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) IA, N
      ELSE IF (IB.LT.MAX(1,N)) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99997) IB, N
      ELSE
C
C        Copy upper triangles of A and B into lower triangles
C        and save the diagonal of B in DL.
         DO 40 J = 1, N
            DO 20 I = 1, J
               A(J,I) = A(I,J)
               B(J,I) = B(I,J)
   20       CONTINUE
            DL(J) = B(J,J)
   40    CONTINUE
C
C        Compute the Cholesky factor L of B.
         CALL F07FDF('Lower',N,B,IB,INFO)
C
         IF (INFO.EQ.0) THEN
C
C           Form the transpose of the upper triangle of inv(L) * A
C           in the lower triangle of array A.
            DO 80 I = 1, N
               Y = B(I,I)
               CALL DGEMV('N',N-I+1,I-1,-ONE,A(I,1),IA,B(I,1),IB,ONE,
     *                    A(I,I),1)
               DO 60 J = I, N
                  A(J,I) = A(J,I)/Y
   60          CONTINUE
   80       CONTINUE
C
C           Form the lower triangle of inv(L) * A * inv(LT)
C           in the lower triangle of array A.
            IB1 = IB
            DO 100 J = 1, N
               CALL DGEMV('N',N-J+1,J-1,-ONE,B(J,1),IB,A(J,1),IA,ONE,
     *                    A(J,J),1)
               IF (J.EQ.N) IB1 = 1
               CALL DTRSV('L','N','N',N-J+1,B(J,J),IB1,A(J,J),1)
  100       CONTINUE
C
C           Restore the diagonal of B and put the diagonal of L into DL.
            CALL DSWAP(N,DL,1,B,IB+1)
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99996)
C           Restore the diagonal of B.
            CALL DCOPY(N,DL,1,B,IB+1)
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99997 FORMAT (1X,'** On entry, IB.lt.max(1,N): IB =',I16,', N =',I16)
99996 FORMAT (1X,'** Matrix B is not positive-definite.')
      END
