      SUBROUTINE F02GAF(JOB,N,A,LDA,W,Z,LDZ,RWORK,WORK,LWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1643 (JUN 1995).
C
C     F02GAF computes all the eigenvalues, the Schur form, or the
C     complete Schur factorization, of a complex general matrix A.
C
C     F02GAF is a driver routine which calls computational routines
C     from LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02GAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDZ, LWORK, N
      CHARACTER         JOB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), W(*), WORK(LWORK), Z(LDZ,*)
      DOUBLE PRECISION  RWORK(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM, TEMP
      DOUBLE PRECISION  AMAX, RMAX, SIGMA
      INTEGER           IBAL, IERR, IHI, ILO, IMAX, INFO, ITAU, IWRK, J,
     *                  LWRK, NREC, I
      LOGICAL           SCALE, WANTS, WANTZ
      CHARACTER         COMPZ, JOB2
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           IDAMAX, IZAMAX, P01ABF
      EXTERNAL          X02AJF, X02AMF, IDAMAX, IZAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02EAZ, F06TFF, P01ABW, P01ABY, ZDSCAL, ZGEBAK,
     *                  ZGEBAL, ZGEHRD, ZHSEQR, ZUNGHR
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
      WANTS = WANTZ .OR. (JOB.EQ.'S') .OR. (JOB.EQ.'s')
C
      IERR = 0
      NREC = 0
      IF ( .NOT. WANTS .AND. JOB.NE.'N' .AND. JOB.NE.'n')
     *    CALL P01ABW(JOB,'JOB',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.MAX(1,N)) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (WANTZ .AND. LDZ.LT.MAX(1,N) .OR. LDZ.LT.1) CALL P01ABY(LDZ,
     *    'LDZ',IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,2*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 80
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
C     Permute the matrix to make it more nearly triangular
C
      CALL ZGEBAL('Permute',N,A,LDA,ILO,IHI,RWORK(IBAL),IERR)
C
C     Reduce matrix to upper Hessenberg form
C
      CALL ZGEHRD(N,ILO,IHI,A,LDA,WORK(ITAU),WORK(IWRK),LWRK,INFO)
C
      IF (WANTZ) THEN
C
C        Copy Householder vectors to Z
C
         CALL F06TFF('Lower',N,N,A,LDA,Z,LDZ)
C
C        Generate unitary matrix Q used in the reduction
C
         CALL ZUNGHR(N,ILO,IHI,Z,LDZ,WORK(ITAU),WORK(IWRK),LWRK,INFO)
      END IF
C
C     Compute eigenvalues, Schur form and Schur vectors, as required
C
      IF (WANTS) THEN
         JOB2 = 'S'
      ELSE
         JOB2 = 'E'
      END IF
      IF (WANTZ) THEN
         COMPZ = 'V'
      ELSE
         COMPZ = 'N'
      END IF
      CALL ZHSEQR(JOB2,COMPZ,N,ILO,IHI,A,LDA,W,Z,LDZ,WORK(IWRK),LWRK,
     *            IERR)
C
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR - 1
         IERR = 2
         NREC = 2
         GO TO 80
      END IF
C
      IF (SCALE) THEN
C
C        Rescale eigenvalues
C
         CALL ZDSCAL(N,ONE/SIGMA,W,1)
         IF (WANTS) THEN
C
C           Rescale Schur form
C
            DO 60 J = 1, N
               CALL ZDSCAL(J,ONE/SIGMA,A(1,J),1)
   60       CONTINUE
         END IF
      END IF
C
      IF (WANTZ) THEN
C
C        Permute Schur vectors to those of the original matrix
C
         CALL ZGEBAK('Permute','Right vectors',N,ILO,IHI,RWORK(IBAL),N,
     *               Z,LDZ,INFO)
C
C        Normalize Schur vectors so that element of largest absolute
C        value is real and positive.
C
         DO 70 J = 1, N
            DO 65 I = 1, N
               RWORK(I) = ABS(Z(I,J))
   65       CONTINUE
            I = IDAMAX(N,RWORK,1)
            IF (DIMAG(Z(I,J)).NE.ZERO .OR. DBLE(Z(I,J)).LT.ZERO) THEN
               TEMP = Z(I,J)/RWORK(I)
               CALL ZSCAL(N,DCONJG(TEMP),Z(1,J),1)
               Z(I,J) = RWORK(I)
               CALL ZSCAL(J-1,DCONJG(TEMP),A(1,J),1)
               IF (J.LT.N) CALL ZSCAL(N-J,TEMP,A(J,J+1),LDA)
            END IF
   70    CONTINUE
      END IF
C
   80 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The QR algorithm failed to converge:',/' ** ',I6,
     *       ' eigenvalues have converged')
      END
