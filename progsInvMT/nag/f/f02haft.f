      SUBROUTINE F02HAF(JOB,UPLO,N,A,LDA,W,RWORK,WORK,LWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1644 (JUN 1995).
C
C     F02HAF computes all the eigenvalues, and optionally all the
C     eigenvectors, of a complex Hermitian matrix A.
C
C     F02HAF is a driver routine which calls computational routines from
C     LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02HAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LWORK, N
      CHARACTER         JOB, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), WORK(LWORK)
      DOUBLE PRECISION  RWORK(*), W(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM, TEMP
      DOUBLE PRECISION  AMAX, RMAX, SIGMA
      INTEGER           I, IE, IERR, IMAX, INFO, IRWK, ITAU, IWRK, J,
     *                  LWRK, NREC
      LOGICAL           SCALE, UPPER, WANTZ
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           IDAMAX, IZAMAX, P01ABF
      EXTERNAL          X02AJF, X02AMF, IDAMAX, IZAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSTERF, F02EAZ, P01ABW, P01ABY, ZDSCAL,
     *                  ZHETRD, ZSTEQR, ZUNGTR
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCONJG, DIMAG, MAX, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(CDUM) = ABS(DBLE(CDUM)) + ABS(DIMAG(CDUM))
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
      IF (LWORK.LT.MAX(1,2*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 100
      END IF
C
C     Check that the diagonal elements of A are real
C
      DO 20 I = 1, N
         IF (DIMAG(A(I,I)).NE.ZERO) THEN
            WRITE (REC,FMT=99997) I
            IERR = 3
            NREC = 1
            GO TO 100
         END IF
   20 CONTINUE
C
C     Get machine constants
C
      RMAX = SQRT(X02AJF()/X02AMF())
C
C     Scale matrix so that maximum magnitude of an element lies
C     in the range ONE to RMAX (magnitude is defined by CABS1)
C
      DO 40 J = 1, N
         IF (UPPER) THEN
            IMAX = IZAMAX(J,A(1,J),1)
         ELSE
            IMAX = J - 1 + IZAMAX(N-J+1,A(J,J),1)
         END IF
         RWORK(J) = CABS1(A(IMAX,J))
   40 CONTINUE
      IMAX = IDAMAX(N,RWORK,1)
      AMAX = RWORK(IMAX)
      CALL F02EAZ(AMAX,ONE,RMAX,SIGMA,SCALE)
      IF (SCALE) THEN
         IF (UPPER) THEN
            DO 60 J = 1, N
               CALL ZDSCAL(J,SIGMA,A(1,J),1)
   60       CONTINUE
         ELSE
            DO 80 J = 1, N
               CALL ZDSCAL(N-J+1,SIGMA,A(J,J),1)
   80       CONTINUE
         END IF
      END IF
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
         GO TO 100
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
C        value is real and positive
C
         DO 90 J = 1, N
            DO 85 I = 1, N
               RWORK(I) = ABS(A(I,J))
   85       CONTINUE
            I = IDAMAX(N,RWORK,1)
            IF (DIMAG(A(I,J)).NE.ZERO .OR. DBLE(A(I,J)).LT.ZERO) THEN
               TEMP = A(I,J)/RWORK(I)
               CALL ZSCAL(N,DCONJG(TEMP),A(1,J),1)
               A(I,J) = RWORK(I)
            END IF
   90    CONTINUE
      END IF
C
  100 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The QR algorithm failed to converge:',/' ** ',I6,
     *' off-diagonal elements of the tridiagonal form have not converged
     *')
99997 FORMAT (' ** A(i,i) has a non-zero imaginary part when i =',I6)
      END
