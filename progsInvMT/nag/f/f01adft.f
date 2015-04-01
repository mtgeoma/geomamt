      SUBROUTINE F01ADF(N,A,IA,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     CHOLINVERSION1
C     The upper triangle of a positive definite symmetric matrix,
C     A, is stored in the upper triangle of an (N+1)*N array
C     A(I,J), I=1,N+1, J=1,N. The Cholesky decomposition A=LU,
C     where U is the transpose of L , is performed and L is stored
C     in the remainder of the array A. The reciprocals of the
C     diagonal elements are stored instead of the elements them-
C     selves. L is then replaced by its inverse and this in turn
C     is replaced by the lower triangle of the inverse of A. A is
C     retained so that the inverse can be subsequently improved.
C     The subroutine will fail if A, modified by the rounding
C     errors, is not positive definite.
C     1st December 1971
C
C     Rewritten to call LAPACK routines SPOTRF/F07FDF and
C     SPOTRI/F07FJF; new IFAIL exit for illegal input
C     parameters inserted; error messages inserted (February 1991).
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01ADF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*)
C     .. Local Scalars ..
      INTEGER           I, IERR, INFO, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F07FDF, F07FJF
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.N+1) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) IA, N
      ELSE
C
C        Copy the upper triangle of A into the lower triangle.
         DO 40 J = 1, N
            DO 20 I = 1, J
               A(J+1,I) = A(I,J)
   20       CONTINUE
   40    CONTINUE
C
C        Factorise A.
         CALL F07FDF('Lower',N,A(2,1),IA,INFO)
         IF (INFO.EQ.0) THEN
C           Invert A.
            CALL F07FJF('Lower',N,A(2,1),IA,INFO)
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99997)
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N = ',I16)
99998 FORMAT (1X,'** On entry, IA.lt.N + 1: IA, N = ',I16,',',I16)
99997 FORMAT (1X,'** Matrix A is not positive-definite.')
      END
