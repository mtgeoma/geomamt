      SUBROUTINE F02GCF(CRIT,N,A,LDA,WL,WU,MEST,M,W,V,LDV,WORK,LWORK,
     *                  RWORK,IWORK,BWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     F02GCF computes selected eigenvalues and eigenvectors of a
C     complex general matrix A.
C
C     F02GCF is a driver routine which calls computational routines
C     from LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02GCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  WL, WU
      INTEGER           IFAIL, LDA, LDV, LWORK, M, MEST, N
      CHARACTER         CRIT
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), V(LDV,MEST), W(*), WORK(LWORK)
      DOUBLE PRECISION  RWORK(*)
      INTEGER           IWORK(*)
      LOGICAL           BWORK(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      DOUBLE PRECISION  ABSW, AMAX, RMAX, SIGMA, WLT, WUT
      INTEGER           I, IBAL, IERR, IHES, IHI, ILO, INFO, IRWRK,
     *                  ITAU, IWRK, J, JJ, LWRK, NREC
      LOGICAL           MODCRI, SCALE
C     .. Local Arrays ..
      COMPLEX*16        VV(1,1)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06UAF, X02AJF, X02AMF
      INTEGER           IDAMAX, P01ABF
      EXTERNAL          DNRM2, F06UAF, X02AJF, X02AMF, IDAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02EAZ, F06TFF, P01ABW, P01ABX, P01ABY, ZDSCAL,
     *                  ZGEBAK, ZGEBAL, ZGEHRD, ZHSEIN, ZHSEQR, ZSCAL,
     *                  ZUNMHR
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCONJG, MAX, MIN, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters
C
      MODCRI = (CRIT.EQ.'M') .OR. (CRIT.EQ.'m')
C
      IERR = 0
      M = 0
      NREC = 0
      IF ( .NOT. MODCRI .AND. CRIT.NE.'R' .AND. CRIT.NE.'r')
     *    CALL P01ABW(CRIT,'CRIT',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.MAX(1,N)) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (WU.LE.WL) CALL P01ABX(WU,'WU',IFAIL,IERR,SRNAME)
      IF (MEST.LT.1) CALL P01ABY(MEST,'MEST',IFAIL,IERR,SRNAME)
      IF (LDV.LT.MAX(1,N)) CALL P01ABY(LDV,'LDV',IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,N*(N+2))) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *    SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 220
      END IF
C
C     Get machine constants
C
      RMAX = SQRT(X02AJF()/X02AMF())
C
C     Scale matrix so that maximum magnitude of an element lies
C     in the range ONE to RMAX
C
      AMAX = F06UAF('Max',N,N,A,LDA,RWORK)
      CALL F02EAZ(AMAX,ONE,RMAX,SIGMA,SCALE)
      IF (SCALE) THEN
         DO 20 J = 1, N
            CALL ZDSCAL(N,SIGMA,A(1,J),1)
   20    CONTINUE
      END IF
C
      IBAL = 1
      IRWRK = IBAL + N
      ITAU = 1
      IWRK = N + ITAU
      IHES = N + IWRK
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
C     Copy Hessenberg matrix to WORK
C
      CALL F06TFF('General',N,N,A,LDA,WORK(IHES),N)
C
C     Compute eigenvalues
C
      CALL ZHSEQR('Eigenvalues only','No Schur vectors',N,ILO,IHI,
     *            WORK(IHES),N,W,VV,1,WORK(IWRK),N,IERR)
C
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR - 1
         IERR = 2
         NREC = 3
         GO TO 180
      END IF
C
C     Select eigenvalues in chosen range
C
      IF (SCALE) THEN
         WUT = WU*SIGMA
         WLT = WL*SIGMA
      ELSE
         WUT = WU
         WLT = WL
      END IF
      M = 0
      DO 40 I = 1, N
         IF (MODCRI) THEN
            ABSW = ABS(W(I))
         ELSE
            ABSW = DBLE(W(I))
         END IF
         BWORK(I) = ABSW .GE. WLT .AND. ABSW .LE. WUT
         IF (BWORK(I)) M = M + 1
   40 CONTINUE
      IF (M.GT.MEST) THEN
         WRITE (REC,FMT=99997) M, MEST
         IERR = 3
         NREC = 4
         GO TO 160
      END IF
C
C     Compute eigenvectors of Hessenberg matrix
C
      CALL ZHSEIN('Right','QR','No initial vectors',BWORK,N,A,LDA,W,VV,
     *            1,V,LDV,MEST,M,WORK(IWRK),RWORK(IRWRK),IWORK,IWORK,
     *            IERR)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99996) IERR
         IERR = 4
         NREC = 3
C
C        Set eigenvectors which have not converged to zero
C
         DO 80 J = 1, M
            IF (IWORK(J).GT.0) THEN
               DO 60 I = 1, N
                  V(I,J) = ZERO
   60          CONTINUE
            END IF
   80    CONTINUE
      END IF
C
C     Transform eigenvectors back to those of original matrix
C
      CALL ZUNMHR('Left','No transpose',N,M,ILO,IHI,A,LDA,WORK(ITAU),V,
     *            LDV,WORK(IWRK),LWRK,INFO)
      CALL ZGEBAK('Both','Right',N,ILO,IHI,RWORK(IBAL),M,V,LDV,INFO)
C
C     Move selected eigenvalues to leading elements of W
C
      J = 0
      DO 100 I = 1, N
         IF (BWORK(I)) THEN
            J = J + 1
            IF (J.NE.I) THEN
               TEMP = W(J)
               W(J) = W(I)
               W(I) = TEMP
            END IF
         END IF
  100 CONTINUE
C
C     Normalize eigenvectors to have unit 2-norm and so that
C     element of largest absolute value is real and positive.
C
      DO 140 J = 1, M
         DO 120 I = 1, N
            RWORK(I) = ABS(V(I,J))
  120    CONTINUE
         JJ = IDAMAX(N,RWORK,1)
         TEMP = V(JJ,J)/(RWORK(JJ)*DNRM2(N,RWORK,1))
         CALL ZSCAL(N,DCONJG(TEMP),V(1,J),1)
         V(JJ,J) = DBLE(V(JJ,J))
  140 CONTINUE
C
  160 CONTINUE
C
      IF (SCALE) THEN
C
C        Rescale eigenvalues
C
         CALL ZDSCAL(N,ONE/SIGMA,W,1)
      END IF
C
  180 CONTINUE
      IF (SCALE) THEN
C
C        Rescale Hessenberg form
C
         DO 200 J = 1, N
            CALL ZDSCAL(MIN(J+1,N),ONE/SIGMA,A(1,J),1)
  200    CONTINUE
      END IF
C
  220 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The QR algorithm failed to converge:',/' ** only',I6,
     *       ' eigenvalues have been computed;',/' ** no eigenvectors',
     *       ' have been computed.')
99997 FORMAT (' ** There are more than MEST eigenvalues in the specifi',
     *       'ed range.',/' ** M (number of eigenvalues in range) =',I6,
     *       '  MEST =',I6,/' ** No eigenvectors have been computed.',
     *       /' ** Rerun with 2nd dimension of V = MEST >= M.')
99996 FORMAT (' ** Inverse iteration failed to compute all the specifi',
     *       'ed eigenvectors.',/' ** The number of eigenvectors which',
     *       ' failed to converge is',I6,/' ** The corresponding colum',
     *       'ns of V are set to zero.')
      END
