      SUBROUTINE F02EAF(JOB,N,A,LDA,WR,WI,Z,LDZ,WORK,LWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1641 (JUN 1995).
C
C     F02EAF computes all the eigenvalues, the Schur form, or the
C     complete Schur factorization, of a real general matrix A.
C
C     F02EAF is a driver routine which calls computational routines
C     from LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDZ, LWORK, N
      CHARACTER         JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WI(*), WORK(LWORK), WR(*), Z(LDZ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAX, RMAX, SIGMA
      INTEGER           I, IBAL, IERR, IHI, ILO, INFO, ITAU, IWRK, J,
     *                  LWRK, NREC
      LOGICAL           SCALE, WANTS, WANTZ
      CHARACTER         COMPZ, JOB2
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  F06RAF, X02AJF, X02AMF
      INTEGER           IDAMAX, P01ABF
      EXTERNAL          F06RAF, X02AJF, X02AMF, IDAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEBAK, DGEBAL, DGEHRD, DHSEQR, DORGHR, DSCAL,
     *                  DSWAP, F02EAZ, F06QFF, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
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
      IF (LWORK.LT.MAX(1,3*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
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
C     Permute the matrix to make it more nearly triangular
C
      CALL DGEBAL('Permute',N,A,LDA,ILO,IHI,WORK(IBAL),INFO)
C
C     Reduce matrix to upper Hessenberg form
C
      CALL DGEHRD(N,ILO,IHI,A,LDA,WORK(ITAU),WORK(IWRK),LWRK,INFO)
C
      IF (WANTZ) THEN
C
C        Copy Householder vectors to Z
C
         CALL F06QFF('Lower',N,N,A,LDA,Z,LDZ)
C
C        Generate orthogonal matrix Q used in the reduction
C
         CALL DORGHR(N,ILO,IHI,Z,LDZ,WORK(ITAU),WORK(IWRK),LWRK,INFO)
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
      CALL DHSEQR(JOB2,COMPZ,N,ILO,IHI,A,LDA,WR,WI,Z,LDZ,WORK(IWRK),
     *            LWRK,IERR)
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
C        Rescale real parts of eigenvalues
C
         CALL DSCAL(N,ONE/SIGMA,WR,1)
C
         IF ( .NOT. WANTS) THEN
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
C                    form, and columns I and I+1 of Z
C
                     CALL DSWAP(I-1,A(1,I),1,A(1,I+1),1)
                     IF (I.LT.N-1) CALL DSWAP(N-I-1,A(I,I+2),LDA,
     *                                        A(I+1,I+2),LDA)
                     A(I,I+1) = A(I+1,I)
                     A(I+1,I) = ZERO
                     IF (WANTZ) CALL DSWAP(N,Z(1,I),1,Z(1,I+1),1)
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
      IF (WANTZ) THEN
C
C        Permute Schur vectors to those of the original matrix
C
         CALL DGEBAK('Permute','Right vectors',N,ILO,IHI,WORK(IBAL),N,Z,
     *               LDZ,INFO)
C
C        Normalize Schur vectors so that element of largest absolute
C        value is positive.
C
         DO 110 J = 1, N
            I = IDAMAX(N,Z(1,J),1)
            IF (Z(I,J).LT.ZERO) THEN
               CALL DSCAL(N,-ONE,Z(1,J),1)
               IF (J.GT.1) THEN
                  CALL DSCAL(J-1,-ONE,A(1,J),1)
                  A(J,J-1) = -A(J,J-1)
               END IF
               IF (J.LT.N) THEN
                  CALL DSCAL(N-J,-ONE,A(J,J+1),LDA)
                  A(J+1,J) = -A(J+1,J)
               END IF
            END IF
  110    CONTINUE
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
