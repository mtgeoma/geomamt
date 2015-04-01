      SUBROUTINE F04ABF(A,IA,B,IB,N,M,C,IC,WKSPCE,BB,IBB,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Accurate solution of a set of real symmetric positive
C     definite linear equations with multiple right hand sides.
C     1st August 1971
C
C     Rewritten to call LAPACK routine SPOTRF/F07FDF;
C     new IFAIL exit inserted for illegal input parameters;
C     error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04ABF')
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IBB, IC, IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), B(IB,*), BB(IBB,*), C(IC,*), WKSPCE(*)
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
      ELSE IF (M.LT.0) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99998) M
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99997) IA, N
      ELSE IF (IB.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99996) IB, N
      ELSE IF (IC.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99995) IC, N
      ELSE IF (IBB.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99994) IBB, N
      ELSE
         IF (N.GT.0) THEN
C           Copy the upper triangle to the lower triangle, and
C           save a copy of the diagonal elements of A.
            DO 40 I = 1, N
               DO 20 J = I + 1, N
                  A(J,I) = A(I,J)
   20          CONTINUE
               WKSPCE(I) = A(I,I)
   40       CONTINUE
C
C           Find the Cholesky factor L.
            CALL F07FDF('Lower',N,A,IA,INFO)
C
            IF (INFO.EQ.0) THEN
C
C              Copy the diagonal elements of A back from WKSPCE, and
C              replace WKSPCE by the reciprocal diagonal elements of L.
   60          DO 80 I = 1, N
                  T = WKSPCE(I)
                  WKSPCE(I) = ONE/A(I,I)
                  A(I,I) = T
   80          CONTINUE
C
               IF (M.GT.0) THEN
C                 Perform iterative refinement on the solution.
                  EPS = X02AJF()
                  IERR2 = 1
                  CALL F04AFF(N,M,A,IA,WKSPCE,B,IB,EPS,C,IC,BB,IBB,NIT,
     *                        IERR2)
                  IF (IERR2.EQ.1) THEN
                     IERR = 2
                     NREC = 2
                     WRITE (P01REC,FMT=99993)
                  END IF
               END IF
C
            ELSE
               IERR = 1
               NREC = 1
               WRITE (P01REC,FMT=99992)
C              Copy the diagonal elements of A back from WKSPCE.
               DO 100 I = 1, N
                  A(I,I) = WKSPCE(I)
  100          CONTINUE
            END IF
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, M.lt.0: M =',I16)
99997 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99996 FORMAT (1X,'** On entry, IB.lt.max(1,N): IB =',I16,', N =',I16)
99995 FORMAT (1X,'** On entry, IC.lt.max(1,N): IC =',I16,', N =',I16)
99994 FORMAT (1X,'** On entry, IBB.lt.max(1,N): IBB =',I16,', N =',I16)
99993 FORMAT (1X,'** Matrix A is too ill-conditioned;',/4X,'iterative ',
     *       'refinement fails to improve the solution.')
99992 FORMAT (1X,'** Matrix A is not positive-definite.')
      END
