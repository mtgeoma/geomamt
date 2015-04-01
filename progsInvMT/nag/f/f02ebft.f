      SUBROUTINE F02EBF(JOB,N,A,LDA,WR,WI,VR,LDVR,VI,LDVI,WORK,LWORK,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     F02EBF computes all the eigenvalues, and optionally all the
C     eigenvectors, of a real general matrix A.
C
C     F02EBF is a driver routine which calls computational routines
C     from LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02EBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDVI, LDVR, LWORK, N
      CHARACTER         JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), VI(LDVI,*), VR(LDVR,*), WI(*),
     *                  WORK(LWORK), WR(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAX, CS, R, RMAX, SCL, SIGMA, SN
      INTEGER           I, IBAL, IERR, IHI, ILO, IMAX, INFO, ITAU, IWRK,
     *                  J, LWRK, M, NREC
      LOGICAL           SCALE, WANTV
      CHARACTER         COMPZ, JOB2
C     .. Local Arrays ..
      LOGICAL           SELECT(1)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BNF, F06RAF, X02AJF, X02AMF
      INTEGER           IDAMAX, P01ABF
      EXTERNAL          DNRM2, F06BNF, F06RAF, X02AJF, X02AMF, IDAMAX,
     *                  P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEBAK, DGEBAL, DGEHRD, DHSEQR, DORGHR,
     *                  DROT, DSCAL, DSWAP, DTREVC, F02EAZ, F06QFF,
     *                  F08HEW, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
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
      IF (WANTV .AND. LDVR.LT.MAX(1,N) .OR. LDVR.LT.1) CALL P01ABY(LDVR,
     *    'LDVR',IFAIL,IERR,SRNAME)
      IF (WANTV .AND. LDVI.LT.MAX(1,N) .OR. LDVI.LT.1) CALL P01ABY(LDVI,
     *    'LDVI',IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,4*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 200
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
      IF (WANTV) THEN
C
C        Copy Householder vectors to VR
C
         CALL F06QFF('Lower',N,N,A,LDA,VR,LDVR)
C
C        Generate the orthogonal matrix Q in VR
C
         CALL DORGHR(N,ILO,IHI,VR,LDVR,WORK(ITAU),WORK(IWRK),LWRK,INFO)
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
      CALL DHSEQR(JOB2,COMPZ,N,ILO,IHI,A,LDA,WR,WI,VR,LDVR,WORK(IWRK),
     *            LWRK,IERR)
C
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR - 1
         IERR = 2
         NREC = 2
         GO TO 200
      END IF
C
      IF (SCALE) THEN
C
C        Rescale real parts of eigenvalues
C
         CALL DSCAL(N,ONE/SIGMA,WR,1)
C
         IF ( .NOT. WANTV) THEN
C
C           Rescale imaginary parts of eigenvalues
C
            CALL DSCAL(N,ONE/SIGMA,WI,1)
         ELSE
C
C           Rescale Schur form
C
            DO 40 J = 1, N
               CALL DSCAL(MIN(J+1,N),ONE/SIGMA,A(1,J),1)
   40       CONTINUE
            IF (SIGMA.GT.ONE) THEN
C
C              Adjust Schur form if rescaling toward underflow has
C              resulted in a 2-by-2 block with A(I+1,I).ne.0 and
C              A(I,I+1).eq.0.
C
               DO 60 I = ILO, IHI - 1
                  IF (A(I+1,I).NE.ZERO .AND. A(I,I+1).EQ.ZERO) THEN
C
C                    Interchange rows and columns I and I+1 in Schur
C                    form, and columns I and I+1 of Schur vectors
C
                     CALL DSWAP(I-1,A(1,I),1,A(1,I+1),1)
                     IF (I.LT.N-1) CALL DSWAP(N-I-1,A(I,I+2),LDA,
     *                                        A(I+1,I+2),LDA)
                     A(I,I+1) = A(I+1,I)
                     A(I+1,I) = ZERO
                     CALL DSWAP(N,VR(1,I),1,VR(1,I+1),1)
                  END IF
   60          CONTINUE
            END IF
C
C           Recompute imaginary parts of eigenvalues after rescaling
C
            DO 80 I = ILO, IHI
               WI(I) = ZERO
   80       CONTINUE
            DO 100 I = ILO, IHI - 1
               IF (A(I+1,I).NE.ZERO) THEN
                  WI(I) = SQRT(ABS(A(I+1,I)))*SQRT(ABS(A(I,I+1)))
                  WI(I+1) = -WI(I)
               END IF
  100       CONTINUE
         END IF
      END IF
C
      IF (WANTV) THEN
C
C        Compute eigenvectors of balanced matrix
C
         IWRK = IBAL + N
         CALL DTREVC('Right','Overwrite',SELECT,N,A,LDA,VI,LDVI,VR,LDVR,
     *               N,M,WORK(IWRK),INFO)
C
C        Transform eigenvectors back to those of original matrix
C
         CALL DGEBAK('Both','Right',N,ILO,IHI,WORK(IBAL),N,VR,LDVR,INFO)
C
C        Normalize eigenvectors to have unit 2-norm and so that
C        element of largest absolute value is real and positive;
C        store real parts in VR and imaginary parts in VI
C
         DO 180 J = 1, N
            IF (WI(J).EQ.ZERO) THEN
               I = IDAMAX(N,VR(1,J),1)
               SCL = DNRM2(N,VR(1,J),1)
               IF (VR(I,J).LT.ZERO)
     *            SCL = -SCL
               CALL DSCAL(N,ONE/SCL,VR(1,J),1)
               DO 120 I = 1, N
                  VI(I,J) = ZERO
  120          CONTINUE
            ELSE IF (WI(J).GT.ZERO) THEN
               SCL = F06BNF(DNRM2(N,VR(1,J),1),DNRM2(N,VR(1,J+1),1))
               CALL DSCAL(N,ONE/SCL,VR(1,J),1)
               CALL DSCAL(N,ONE/SCL,VR(1,J+1),1)
               DO 140 I = 1, N
                  WORK(I) = VR(I,J)**2 + VR(I,J+1)**2
  140          CONTINUE
               IMAX = IDAMAX(N,WORK,1)
               CALL F08HEW(VR(IMAX,J),VR(IMAX,J+1),CS,SN,R)
               IF (R.LT.ZERO) THEN
                  CS = -CS
                  SN = -SN
               END IF
               CALL DROT(N,VR(1,J),1,VR(1,J+1),1,CS,SN)
               VR(IMAX,J+1) = ZERO
               CALL DCOPY(N,VR(1,J+1),1,VI(1,J),1)
               CALL DCOPY(N,VR(1,J),1,VR(1,J+1),1)
               DO 160 I = 1, N
                  VI(I,J+1) = -VI(I,J)
  160          CONTINUE
            END IF
  180    CONTINUE
      END IF
C
  200 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The QR algorithm failed to converge:',/' ** ',I6,
     *       ' eigenvalues have converged')
      END
