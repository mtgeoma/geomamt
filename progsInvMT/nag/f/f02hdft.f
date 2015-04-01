      SUBROUTINE F02HDF(ITYPE,JOB,UPLO,N,A,LDA,B,LDB,W,RWORK,WORK,LWORK,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     F02HDF computes all the eigenvalues, and optionally all the
C     eigenvectors, of a complex Hermitian-definite generalized
C     eigenproblem, of one of the types:
C
C       1. A*z = lambda*B*z
C       2. A*B*z = lambda*z
C       3. B*A*z = lambda*z.
C
C     Here A and B are Hermitian and B is also positive-definite.
C
C     F02HDF is a driver routine which calls computational routines from
C     LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02HDF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, ITYPE, LDA, LDB, LWORK, N
      CHARACTER         JOB, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), WORK(LWORK)
      DOUBLE PRECISION  RWORK(*), W(*)
C     .. Local Scalars ..
      INTEGER           I, IE, IERR, INFO, IRWK, ITAU, IWRK, LWRK, NREC
      LOGICAL           UPPER, WANTZ
      CHARACTER         TRANS
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DSTERF, P01ABW, P01ABY, ZHEGST, ZHETRD, ZPOTRF,
     *                  ZSTEQR, ZTRMM, ZTRSM, ZUNGTR
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, MAX
C     .. Executable Statements ..
C
C     Test the input parameters
C
      WANTZ = (JOB.EQ.'V') .OR. (JOB.EQ.'v')
      UPPER = (UPLO.EQ.'U') .OR. (UPLO.EQ.'u')
C
      IERR = 0
      NREC = 0
      IF (ITYPE.LT.1 .OR. ITYPE.GT.3) CALL P01ABY(ITYPE,'ITYPE',IFAIL,
     *    IERR,SRNAME)
      IF ( .NOT. WANTZ .AND. JOB.NE.'N' .AND. JOB.NE.'n')
     *    CALL P01ABW(JOB,'JOB',IFAIL,IERR,SRNAME)
      IF ( .NOT. UPPER .AND. UPLO.NE.'L' .AND. UPLO.NE.'l')
     *    CALL P01ABW(UPLO,'UPLO',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.MAX(1,N)) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (LDB.LT.MAX(1,N)) CALL P01ABY(LDB,'LDB',IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,2*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 60
      END IF
C
C     Check that the diagonal elements of A are real
C
      DO 20 I = 1, N
         IF (DIMAG(A(I,I)).NE.ZERO) THEN
            WRITE (REC,FMT=99996) I
            IERR = 4
            NREC = 1
            GO TO 60
         END IF
   20 CONTINUE
C
C     Check that the diagonal elements of B are real
C
      DO 40 I = 1, N
         IF (DIMAG(B(I,I)).NE.ZERO) THEN
            WRITE (REC,FMT=99995) I
            IERR = 5
            NREC = 1
            GO TO 60
         END IF
   40 CONTINUE
C
C     Form a Cholesky factorization of B.
C
      CALL ZPOTRF(UPLO,N,B,LDB,IERR)
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99997) IERR
         IERR = 3
         NREC = 2
         GO TO 60
      END IF
C
C     Transform problem to standard eigenvalue problem
C
      CALL ZHEGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
C
      IE = 1
      IRWK = IE + N
      ITAU = 1
      IWRK = ITAU + N
      LWRK = LWORK - IWRK + 1
C
C     Reduce Hermitian matrix to real tridiagonal form
C
      CALL ZHETRD(UPLO,N,A,LDA,W,RWORK(IE),WORK(ITAU),WORK(IWRK),LWRK,
     *            INFO)
C
      IF ( .NOT. WANTZ) THEN
C
C        Compute eigenvalues only
C
         CALL DSTERF(N,W,RWORK(IE),IERR)
      ELSE
C
C        Generate unitary matrix Q used in the reduction
C
         CALL ZUNGTR(UPLO,N,A,LDA,WORK(ITAU),WORK(IWRK),LWRK,INFO)
C
C        Compute eigenvalues and eigenvectors
C
         CALL ZSTEQR(JOB,N,W,RWORK(IE),A,LDA,RWORK(IRWK),IERR)
      END IF
C
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR
         IERR = 2
         NREC = 2
         GO TO 60
      END IF
C
      IF (WANTZ) THEN
C
C        Transform eigenvectors back to those of the original problem
C
         IF (ITYPE.EQ.1 .OR. ITYPE.EQ.2) THEN
            IF (UPPER) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'C'
            END IF
            CALL ZTRSM('Left',UPLO,TRANS,'Non-unit',N,N,ONE,B,LDB,A,LDA)
         ELSE IF (ITYPE.EQ.3) THEN
            IF (UPPER) THEN
               TRANS = 'C'
            ELSE
               TRANS = 'N'
            END IF
            CALL ZTRMM('Left',UPLO,TRANS,'Non-unit',N,N,ONE,B,LDB,A,LDA)
         END IF
      END IF
C
   60 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The QR algorithm failed to converge:',/' ** ',I6,
     *' off-diagonal elements of the tridiagonal form have not converged
     *')
99997 FORMAT (' ** B is not positive-definite:',/' ** Its leading mino',
     *       'r of order ',I6,' is not positive-definite.')
99996 FORMAT (' ** A(i,i) has a non-zero imaginary part when i =',I6)
99995 FORMAT (' ** B(i,i) has a non-zero imaginary part when i =',I6)
      END
