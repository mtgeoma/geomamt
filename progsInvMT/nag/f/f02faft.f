      SUBROUTINE F02FAF(JOB,UPLO,N,A,LDA,W,WORK,LWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1642 (JUN 1995).
C
C     F02FAF computes all the eigenvalues, and optionally all the
C     eigenvectors, of a real symmetric matrix A.
C
C     F02FAF is a driver routine which calls computational routines from
C     LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02FAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LWORK, N
      CHARACTER         JOB, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), W(*), WORK(LWORK)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAX, RMAX, SIGMA
      INTEGER           I, IE, IERR, INFO, ITAU, IWRK, J, LWRK, NREC
      LOGICAL           SCALE, UPPER, WANTZ
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F06RCF, X02AJF, X02AMF
      INTEGER           IDAMAX, P01ABF
      EXTERNAL          F06RCF, X02AJF, X02AMF, IDAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DORGTR, DSCAL, DSTEQR, DSTERF, DSYTRD, F02EAZ,
     *                  P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters
C
      WANTZ = (JOB.EQ.'V') .OR. (JOB.EQ.'v')
      UPPER = (UPLO.EQ.'U') .OR. (UPLO.EQ.'u')
C
      IERR = 0
      NREC = 0
      IF ( .NOT. WANTZ .AND. JOB.NE.'N' .AND. JOB.NE.'n')
     *    CALL P01ABW(JOB,'JOB',IFAIL,IERR,SRNAME)
      IF ( .NOT. UPPER .AND. UPLO.NE.'L' .AND. UPLO.NE.'l')
     *    CALL P01ABW(UPLO,'UPLO',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.MAX(1,N)) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,3*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 60
      END IF
C
C     Get machine constants
C
      RMAX = SQRT(X02AJF()/X02AMF())
C
C     Scale matrix so that maximum magnitude of an element lies
C     in the range ONE to RMAX
C
      AMAX = F06RCF('Max',UPLO,N,A,LDA,WORK)
      CALL F02EAZ(AMAX,ONE,RMAX,SIGMA,SCALE)
      IF (SCALE) THEN
         IF (UPPER) THEN
            DO 20 J = 1, N
               CALL DSCAL(J,SIGMA,A(1,J),1)
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               CALL DSCAL(N-J+1,SIGMA,A(J,J),1)
   40       CONTINUE
         END IF
      END IF
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
         GO TO 60
      END IF
C
      IF (SCALE) THEN
C
C        Rescale eigenvalues
C
         CALL DSCAL(N,ONE/SIGMA,W,1)
      END IF
C
      IF (WANTZ) THEN
C
C        Normalize eigenvectors so that element of largest absolute
C        value is positive.
C
         DO 50 J = 1, N
            I = IDAMAX(N,A(1,J),1)
            IF (A(I,J).LT.ZERO) CALL DSCAL(N,-ONE,A(1,J),1)
   50    CONTINUE
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
      END
