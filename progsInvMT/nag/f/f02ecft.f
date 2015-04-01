      SUBROUTINE F02ECF(CRIT,N,A,LDA,WL,WU,MEST,M,WR,WI,VR,LDVR,VI,LDVI,
     *                  WORK,LWORK,IWORK,BWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     F02ECF computes selected eigenvalues and eigenvectors of a real
C     general matrix A.
C
C     F02ECF is a driver routine which calls computational routines
C     from LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02ECF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  WL, WU
      INTEGER           IFAIL, LDA, LDVI, LDVR, LWORK, M, MEST, N
      CHARACTER         CRIT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), VI(LDVI,MEST), VR(LDVR,MEST), WI(*),
     *                  WORK(LWORK), WR(*)
      INTEGER           IWORK(*)
      LOGICAL           BWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSW, AMAX, CS, RMAX, SIGMA, SN, TEMP, WLT, WUT
      INTEGER           I, IBAL, IERR, IHES, IHI, ILO, INFO, ITAU, IWRK,
     *                  J, JJ, LWRK, NREC
      LOGICAL           MODCRI, SCALE
C     .. Local Arrays ..
      DOUBLE PRECISION  VV(1,1)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BNF, F06RAF, X02AJF, X02AMF
      INTEGER           IDAMAX, P01ABF
      EXTERNAL          DNRM2, F06BNF, F06RAF, X02AJF, X02AMF, IDAMAX,
     *                  P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEBAK, DGEBAL, DGEHRD, DHSEIN, DHSEQR,
     *                  DORMHR, DROT, DSCAL, F02EAZ, F06QFF, F08HEW,
     *                  P01ABW, P01ABX, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
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
      IF (LDVR.LT.MAX(1,N)) CALL P01ABY(LDVR,'LDVR',IFAIL,IERR,SRNAME)
      IF (LDVI.LT.MAX(1,N)) CALL P01ABY(LDVI,'LDVI',IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,N*(N+4))) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *    SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 260
      END IF
C
C     Get machine constants
C
      RMAX = SQRT(X02AJF()/X02AMF())
C
C     Scale matrix so that maximum magnitude of an element lies
C     in the range ONE to RMAX
C
      AMAX = F06RAF('Max',N,N,A,LDA,WORK)
      CALL F02EAZ(AMAX,ONE,RMAX,SIGMA,SCALE)
      IF (SCALE) THEN
         DO 20 J = 1, N
            CALL DSCAL(N,SIGMA,A(1,J),1)
   20    CONTINUE
      END IF
C
      IBAL = 1
      ITAU = N + IBAL
      IWRK = N + ITAU
      IHES = N + IWRK
      LWRK = LWORK - IWRK + 1
C
C     Balance the matrix
C
      CALL DGEBAL('Both',N,A,LDA,ILO,IHI,WORK(IBAL),INFO)
C
C     Reduce to upper Hessenberg form
C
      CALL DGEHRD(N,ILO,IHI,A,LDA,WORK(ITAU),WORK(IWRK),LWRK,INFO)
C
C     Copy Hessenberg matrix to WORK
C
      CALL F06QFF('General',N,N,A,LDA,WORK(IHES),N)
C
C     Compute eigenvalues
C
      CALL DHSEQR('Eigenvalues only','No Schur vectors',N,ILO,IHI,
     *            WORK(IHES),N,WR,WI,VV,1,WORK(IWRK),N,IERR)
C
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR - 1
         IERR = 2
         NREC = 3
         GO TO 220
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
            IF (WI(I).EQ.ZERO) THEN
               ABSW = ABS(WR(I))
            ELSE
               ABSW = F06BNF(WR(I),WI(I))
            END IF
         ELSE
            ABSW = WR(I)
         END IF
         BWORK(I) = ABSW .GE. WLT .AND. ABSW .LE. WUT
         IF (BWORK(I)) M = M + 1
   40 CONTINUE
      IF (M.GT.MEST) THEN
         WRITE (REC,FMT=99997) M, MEST
         IERR = 3
         NREC = 4
         GO TO 200
      END IF
C
C     Compute eigenvectors of Hessenberg matrix
C
      CALL DHSEIN('Right','QR','No initial vectors',BWORK,N,A,LDA,WR,WI,
     *            VI,LDVI,VR,LDVR,MEST,M,WORK(IWRK),IWORK,IWORK,IERR)
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
                  VR(I,J) = ZERO
                  VI(I,J) = ZERO
   60          CONTINUE
            END IF
   80    CONTINUE
      END IF
C
C     Transform eigenvectors back to those of original matrix
C
      CALL DORMHR('Left','No transpose',N,M,ILO,IHI,A,LDA,WORK(ITAU),VR,
     *            LDVR,WORK(IWRK),LWRK,INFO)
      CALL DGEBAK('Both','Right',N,ILO,IHI,WORK(IBAL),M,VR,LDVR,INFO)
C
C     Move selected eigenvalues to leading elements of WR and WI
C     (DHSEIN has modified BWORK so that BWORK(i) = .FALSE. if it
C     corresponds to the second eigenvalue in a complex conjugate
C     pair)
C
      J = 0
      DO 100 I = 1, N
         IF (BWORK(I)) THEN
            IF (WI(I).GT.ZERO) BWORK(I+1) = .TRUE.
            J = J + 1
            IF (J.NE.I) THEN
               TEMP = WR(J)
               WR(J) = WR(I)
               WR(I) = TEMP
               TEMP = WI(J)
               WI(J) = WI(I)
               WI(I) = TEMP
            END IF
         END IF
  100 CONTINUE
C
C     Normalize eigenvectors to have unit 2-norm and so that
C     element of largest absolute value is real and positive;
C     store real parts in VR and imaginary parts in VI
C
      DO 180 J = 1, M
         IF (WI(J).EQ.ZERO) THEN
            I = IDAMAX(N,VR(1,J),1)
            TEMP = DNRM2(N,VR(1,J),1)
            IF (VR(I,J).LT.ZERO) TEMP = -TEMP
            IF (TEMP.NE.ZERO) THEN
               CALL DSCAL(N,ONE/TEMP,VR(1,J),1)
            END IF
            DO 120 I = 1, N
               VI(I,J) = ZERO
  120       CONTINUE
         ELSE IF (WI(J).GT.ZERO) THEN
            TEMP = F06BNF(DNRM2(N,VR(1,J),1),DNRM2(N,VR(1,J+1),1))
            IF (TEMP.GT.ZERO) THEN
               CALL DSCAL(N,ONE/TEMP,VR(1,J),1)
               CALL DSCAL(N,ONE/TEMP,VR(1,J+1),1)
            END IF
            DO 140 I = 1, N
               WORK(I) = VR(I,J)**2 + VR(I,J+1)**2
  140       CONTINUE
            JJ = IDAMAX(N,WORK,1)
            CALL F08HEW(VR(JJ,J),VR(JJ,J+1),CS,SN,TEMP)
            IF (TEMP.LT.ZERO) THEN
               CS = -CS
               SN = -SN
            END IF
            CALL DROT(N,VR(1,J),1,VR(1,J+1),1,CS,SN)
            VR(JJ,J+1) = ZERO
            CALL DCOPY(N,VR(1,J+1),1,VI(1,J),1)
            CALL DCOPY(N,VR(1,J),1,VR(1,J+1),1)
            DO 160 I = 1, N
               VI(I,J+1) = -VI(I,J)
  160       CONTINUE
         END IF
  180 CONTINUE
C
  200 CONTINUE
      IF (SCALE) THEN
C
C        Rescale eigenvalues
C
         CALL DSCAL(N,ONE/SIGMA,WR,1)
         CALL DSCAL(N,ONE/SIGMA,WI,1)
      END IF
C
  220 CONTINUE
      IF (SCALE) THEN
C
C        Rescale Hessenberg form
C
         DO 240 J = 1, N
            CALL DSCAL(MIN(J+1,N),ONE/SIGMA,A(1,J),1)
  240    CONTINUE
      END IF
C
  260 CONTINUE
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
     *       /' ** Rerun with 2nd dimension of VR and VI = MEST >= M.')
99996 FORMAT (' ** Inverse iteration failed to compute all the specifi',
     *       'ed eigenvectors.',/' ** The number of eigenvectors which',
     *       ' failed to converge is',I6,/' ** The corresponding colum',
     *       'ns of VR and VI are set to zero.')
      END
