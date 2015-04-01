      SUBROUTINE F04ASF(A,IA,B,N,C,WK1,WK2,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1004 (JUN 1993).
C
C     Accurate solution of a set of real symmetric positive
C     definite linear equations with one right hand side.
C     1st April 1973
C
C     Rewritten to call LAPACK routine SPOTRF/F07FDF;
C     new IFAIL exit inserted for illegal input parameters;
C     error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04ASF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), B(*), C(*), WK1(*), WK2(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, T
      INTEGER           I, IERR, IERR2, INFO, J, NIT, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04AFF, F07FDF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99998) IA, N
      ELSE
         IF (N.GT.0) THEN
C           Copy the upper triangle to the lower triangle, and
C           save a copy of the diagonal elements of A.
            DO 40 I = 1, N
               DO 20 J = I + 1, N
                  A(J,I) = A(I,J)
   20          CONTINUE
               WK1(I) = A(I,I)
   40       CONTINUE
C
C           Find the Cholesky factor L.
            CALL F07FDF('Lower',N,A,IA,INFO)
C
            IF (INFO.EQ.0) THEN
C
C              Copy the diagonal elements of A back from WK1, and
C              replace WK1 by the reciprocal diagonal elements of L.
               DO 80 I = 1, N
                  T = WK1(I)
                  WK1(I) = ONE/A(I,I)
                  A(I,I) = T
   80          CONTINUE
C
C              Perform iterative refinement on the solution.
               EPS = X02AJF()
               IERR2 = 1
               CALL F04AFF(N,1,A,IA,WK1,B,N,EPS,C,N,WK2,N,NIT,IERR2)
               IF (IERR2.EQ.1) THEN
                  IERR = 2
                  NREC = 2
                  WRITE (P01REC,FMT=99997)
               END IF
C
            ELSE
               IERR = 1
               NREC = 1
               WRITE (P01REC,FMT=99996)
C              Copy the diagonal elements of A back from WK1.
               DO 100 I = 1, N
                  A(I,I) = WK1(I)
  100          CONTINUE
            END IF
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99997 FORMAT (1X,'** Matrix A is too ill-conditioned;',/4X,'iterative ',
     *       'refinement fails to improve the solution.')
99996 FORMAT (1X,'** Matrix A is not positive-definite.')
      END
