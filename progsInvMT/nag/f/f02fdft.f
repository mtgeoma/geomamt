      SUBROUTINE F02FDF(ITYPE,JOB,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     F02FDF computes all the eigenvalues, and optionally all the
C     eigenvectors, of a real symmetric-definite generalized
C     eigenproblem, of one of the types:
C
C       1. A*z = lambda*B*z
C       2. A*B*z = lambda*z
C       3. B*A*z = lambda*z.
C
C     Here A and B are symmetric and B is also positive-definite.
C
C     F02FDF is a driver routine which calls computational routines from
C     LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02FDF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, ITYPE, LDA, LDB, LWORK, N
      CHARACTER         JOB, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), W(*), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           IE, IERR, INFO, ITAU, IWRK, LWRK, NREC
      LOGICAL           UPPER, WANTZ
      CHARACTER         TRANS
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DORGTR, DPOTRF, DSTEQR, DSTERF, DSYGST, DSYTRD,
     *                  DTRMM, DTRSM, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
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
      IF (LWORK.LT.MAX(1,3*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 20
      END IF
C
C     Form a Cholesky factorization of B.
C
      CALL DPOTRF(UPLO,N,B,LDB,IERR)
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99997) IERR
         IERR = 3
         NREC = 2
         GO TO 20
      END IF
C
C     Transform problem to standard eigenvalue problem
C
      CALL DSYGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
C
      IE = 1
      ITAU = IE + N
      IWRK = ITAU + N
      LWRK = LWORK - IWRK + 1
C
C     Reduce symmetric matrix to tridiagonal form
C
      CALL DSYTRD(UPLO,N,A,LDA,W,WORK(IE),WORK(ITAU),WORK(IWRK),LWRK,
     *            INFO)
C
      IF ( .NOT. WANTZ) THEN
C
C        Compute eigenvalues only
C
         CALL DSTERF(N,W,WORK(IE),IERR)
      ELSE
C
C        Generate orthogonal matrix Q used in the reduction
C
         CALL DORGTR(UPLO,N,A,LDA,WORK(ITAU),WORK(IWRK),LWRK,INFO)
C
C        Compute eigenvalues and eigenvectors
C
         CALL DSTEQR(JOB,N,W,WORK(IE),A,LDA,WORK(ITAU),IERR)
      END IF
C
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR
         IERR = 2
         NREC = 2
         GO TO 20
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
               TRANS = 'T'
            END IF
            CALL DTRSM('Left',UPLO,TRANS,'Non-unit',N,N,ONE,B,LDB,A,LDA)
         ELSE IF (ITYPE.EQ.3) THEN
            IF (UPPER) THEN
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            END IF
            CALL DTRMM('Left',UPLO,TRANS,'Non-unit',N,N,ONE,B,LDB,A,LDA)
         END IF
      END IF
C
   20 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The QR algorithm failed to converge:',/' ** ',I6,
     *' off-diagonal elements of the tridiagonal form have not converged
     *')
99997 FORMAT (' ** B is not positive-definite:',/' ** Its leading mino',
     *       'r of order ',I6,' is not positive-definite.')
      END
