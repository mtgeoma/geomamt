      SUBROUTINE F03AEF(N,A,IA,P,D1,ID,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     CHOLDET1
C     The upper triangle of a positive definite symmetric matrix,
C     A, is stored in the upper triangle of an N*N array A(I,J),
C     I=1,N, J=1,N. The Cholesky decomposition A = LU, where
C     U is the transpose of L, is performed and stored in the
C     remainder of the array A except for the reciprocals of the
C     diagonal elements which are stored in P(I), I = 1,N,
C     instead of the elements themselves. A is retained so that
C     the solution obtained can be subsequently improved. The
C     determinant, D1 * 2.0**ID, of A is also computed. The
C     subroutine will fail if A, modified by the rounding errors,
C     is not positive definite.
C     1st December 1971
C
C     Rewritten to call LAPACK routine SPOTRF/F07ADF; new IFAIL
C     exit inserted for illegal input parameters; error messages
C     inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03AEF')
      DOUBLE PRECISION  ONE, ZERO, SIXTEN, SXTNTH
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0,SIXTEN=16.0D+0,
     *                  SXTNTH=6.25D-2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D1
      INTEGER           IA, ID, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), P(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           I, IERR, INFO, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F07FDF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) IA, N
      ELSE
C
C        Copy the upper triangle to the lower triangle, and
C        save a copy of the diagonal elements of A.
         DO 40 I = 1, N
            DO 20 J = I + 1, N
               A(J,I) = A(I,J)
   20       CONTINUE
            P(I) = A(I,I)
   40    CONTINUE
C
         CALL F07FDF('Lower',N,A,IA,INFO)
C
         IF (INFO.EQ.0) THEN
C           Compute the determinant as the product of the squares
C           of the diagonal elements of the Cholesky factor L.
C           Store the determinant as D = D1*2.0**ID.
            D1 = ONE
            ID = 0
            DO 100 I = 1, N
               D1 = D1*A(I,I)*A(I,I)
   60          IF (D1.GE.ONE) THEN
                  D1 = D1*SXTNTH
                  ID = ID + 4
                  GO TO 60
               END IF
   80          IF (D1.LT.SXTNTH) THEN
                  D1 = D1*SIXTEN
                  ID = ID - 4
                  GO TO 80
               END IF
  100       CONTINUE
C
C           Copy the diagonal elements of A back from P, and
C           replace P by the reciprocal diagonal elements of L.
  120       DO 140 I = 1, N
               T = P(I)
               P(I) = ONE/A(I,I)
               A(I,I) = T
  140       CONTINUE
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99997)
            D1 = ZERO
            ID = 0
C           Copy the diagonal elements of A back from P.
            DO 160 I = 1, N
               A(I,I) = P(I)
  160       CONTINUE
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99997 FORMAT (1X,'** Matrix A is not positive-definite.')
      END
