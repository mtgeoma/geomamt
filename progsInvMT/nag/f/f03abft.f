      SUBROUTINE F03ABF(A,IA,N,DET,WKSPCE,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Determinant of real symmetric positive definite matrix
C     1st August 1971
C
C     Rewritten to call LAPACK routine SPOTRF/F07FDF;
C     new IFAIL exit inserted for illegal input parameters;
C     error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03ABF')
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DET
      INTEGER           IA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), WKSPCE(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  OFLOW, T, UFLOW
      INTEGER           I, IERR, INFO, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F07FDF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 4
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 4
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
            WKSPCE(I) = A(I,I)
   40    CONTINUE
C
         CALL F07FDF('Lower',N,A,IA,INFO)
C
         IF (INFO.EQ.0) THEN
C
C           Compute the determinant as the product of the squares
C           of the diagonal elements of the Cholesky factor L.
            UFLOW = X02AMF()
            OFLOW = ONE/UFLOW
            DET = ONE
            DO 60 I = 1, N
               T = A(I,I)
               IF (T.GE.ONE) THEN
                  IF (DET.GT.(OFLOW/T)/T) THEN
                     IERR = 2
                     NREC = 1
                     WRITE (P01REC,FMT=99997)
                     DET = OFLOW
                     GO TO 80
                  ELSE
                     DET = DET*T*T
                  END IF
               ELSE
                  IF (DET.LT.(UFLOW/T)/T) THEN
                     IERR = 3
                     NREC = 1
                     WRITE (P01REC,FMT=99996)
                     DET = ZERO
                     GO TO 80
                  ELSE
                     DET = DET*T*T
                  END IF
               END IF
   60       CONTINUE
C
C           Copy the diagonal elements of A back from WKSPCE, and
C           replace WKSPCE by the reciprocal diagonal elements of L.
   80       DO 100 I = 1, N
               T = WKSPCE(I)
               WKSPCE(I) = ONE/A(I,I)
               A(I,I) = T
  100       CONTINUE
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99995)
            DET = ZERO
C           Copy the diagonal elements of A back from WKSPCE.
            DO 120 I = 1, N
               A(I,I) = WKSPCE(I)
  120       CONTINUE
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99997 FORMAT (1X,'** The value of the determinant is too large to be s',
     *       'tored.')
99996 FORMAT (1X,'** The value of the determinant is too small to be s',
     *       'tored.')
99995 FORMAT (1X,'** Matrix A is not positive-definite.')
      END
