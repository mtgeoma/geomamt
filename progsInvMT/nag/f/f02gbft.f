      SUBROUTINE F02GBF(JOB,N,A,LDA,W,V,LDV,RWORK,WORK,LWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     F02GBF computes all the eigenvalues, and optionally all the
C     eigenvectors, of a complex general matrix A.
C
C     F02GBF is a driver routine which calls computational routines
C     from LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02GBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDV, LWORK, N
      CHARACTER         JOB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), V(LDV,*), W(*), WORK(LWORK)
      DOUBLE PRECISION  RWORK(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM
      DOUBLE PRECISION  AMAX, RMAX, SCL, SIGMA
      INTEGER           I, IBAL, IERR, IHI, ILO, IMAX, INFO, IRWK, ITAU,
     *                  IWRK, J, LWRK, M, NREC
      LOGICAL           SCALE, WANTV
      CHARACTER         COMPZ, JOB2
C     .. Local Arrays ..
      COMPLEX*16        VL(1,1)
      LOGICAL           SELECT(1)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DZNRM2, X02AJF, X02AMF
      INTEGER           IDAMAX, IZAMAX, P01ABF
      EXTERNAL          DZNRM2, X02AJF, X02AMF, IDAMAX, IZAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02EAZ, F06TFF, P01ABW, P01ABY, ZDSCAL, ZGEBAK,
     *                  ZGEBAL, ZGEHRD, ZHSEQR, ZSCAL, ZTREVC, ZUNGHR
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
      WANTV = (JOB.EQ.'V') .OR. (JOB.EQ.'v')
C
      IERR = 0
      NREC = 0
      IF ( .NOT. WANTV .AND. JOB.NE.'N' .AND. JOB.NE.'n')
     *    CALL P01ABW(JOB,'JOB',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.MAX(1,N)) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (WANTV .AND. LDV.LT.MAX(1,N) .OR. LDV.LT.1) CALL P01ABY(LDV,
     *    'LDV',IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,2*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 120
      END IF
C
C     Get machine constants
C
      RMAX = SQRT(X02AJF()/X02AMF())
C
C     Scale matrix so that maximum magnitude of an element lies
C     in the range ONE to RMAX (magnitude is defined by CABS1)
C
      DO 20 J = 1, N
         IMAX = IZAMAX(N,A(1,J),1)
         RWORK(J) = CABS1(A(IMAX,J))
   20 CONTINUE
      IMAX = IDAMAX(N,RWORK,1)
      AMAX = RWORK(IMAX)
      CALL F02EAZ(AMAX,ONE,RMAX,SIGMA,SCALE)
      IF (SCALE) THEN
         DO 40 J = 1, N
            CALL ZDSCAL(N,SIGMA,A(1,J),1)
   40    CONTINUE
      END IF
C
      IBAL = 1
      ITAU = 1
      IWRK = N + ITAU
      LWRK = LWORK - IWRK + 1
C
C     Balance the matrix
C
      CALL ZGEBAL('Both',N,A,LDA,ILO,IHI,RWORK(IBAL),INFO)
C
C     Reduce to upper Hessenberg form
C
      CALL ZGEHRD(N,ILO,IHI,A,LDA,WORK(ITAU),WORK(IWRK),LWRK,INFO)
C
      IF (WANTV) THEN
C
C        Copy Householder vectors to V
C
         CALL F06TFF('Lower',N,N,A,LDA,V,LDV)
C
C        Generate the unitary matrix Q in V
C
         CALL ZUNGHR(N,ILO,IHI,V,LDV,WORK(ITAU),WORK(IWRK),LWRK,INFO)
      END IF
C
C     Compute eigenvalues, and the complete Schur factorization if
C     eigenvectors are required
C
      IF (WANTV) THEN
         JOB2 = 'S'
         COMPZ = 'V'
      ELSE
         JOB2 = 'E'
         COMPZ = 'N'
      END IF
      CALL ZHSEQR(JOB2,COMPZ,N,ILO,IHI,A,LDA,W,V,LDV,WORK(IWRK),LWRK,
     *            IERR)
C
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR - 1
         IERR = 2
         NREC = 2
         GO TO 120
      END IF
C
      IF (SCALE) THEN
C
C        Rescale eigenvalues
C
         CALL ZDSCAL(N,ONE/SIGMA,W,1)
C
         IF (WANTV) THEN
C
C           Rescale Schur form
C
            DO 60 J = 1, N
               CALL ZDSCAL(J,ONE/SIGMA,A(1,J),1)
   60       CONTINUE
         END IF
      END IF
C
      IF (WANTV) THEN
C
C        Compute eigenvectors of balanced matrix
C
         IWRK = 1
         IRWK = IBAL + N
         CALL ZTREVC('Right','Overwrite',SELECT,N,A,LDA,VL,1,V,LDV,N,M,
     *               WORK(IWRK),RWORK(IRWK),INFO)
C
C        Transform eigenvectors back to those of original matrix
C
         CALL ZGEBAK('Both','Right',N,ILO,IHI,RWORK(IBAL),N,V,LDV,INFO)
C
C        Normalize eigenvectors to have unit 2-norm and so that
C        element of largest absolute value is real and positive
C
         DO 100 J = 1, N
            SCL = DZNRM2(N,V(1,J),1)
            DO 80 I = 1, N
               RWORK(I) = ABS(V(I,J))
   80       CONTINUE
            IMAX = IDAMAX(N,RWORK,1)
            SCL = SCL*RWORK(IMAX)
            CALL ZSCAL(N,DCONJG(V(IMAX,J))/SCL,V(1,J),1)
            V(IMAX,J) = DBLE(V(IMAX,J))
  100    CONTINUE
      END IF
C
  120 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The QR algorithm failed to converge:',/' ** ',I6,
     *       ' eigenvalues have converged')
      END
