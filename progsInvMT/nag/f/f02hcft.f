      SUBROUTINE F02HCF(JOB,RANGE,UPLO,N,A,LDA,WL,WU,IL,IU,MEST,M,W,Z,
     *                  LDZ,WORK,LWORK,RWORK,IWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     F02HCF computes selected eigenvalues, and optionally the
C     corresponding eigenvectors, of a complex Hermitian matrix A.
C
C     F02HCF is a driver routine which calls computational routines from
C     LAPACK in Chapter F08.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02HCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  WL, WU
      INTEGER           IFAIL, IL, IU, LDA, LDZ, LWORK, M, MEST, N
      CHARACTER         JOB, RANGE, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), WORK(LWORK), Z(LDZ,MEST)
      DOUBLE PRECISION  RWORK(*), W(*)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      DOUBLE PRECISION  ABSTOL, AMAX, RMAX, SIGMA, TMP1, WLT, WUT
      INTEGER           I, IBL, ID, IE, IERR, IFA, IMAX, INFO, IRWRK,
     *                  ISP, ITAU, IWO, IWRK, J, JJ, LWRK, NREC, NSPLIT
      LOGICAL           SCALE, UPPER, VALEIG, WANTZ
      CHARACTER         ORDER
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           IDAMAX, IZAMAX, P01ABF
      EXTERNAL          X02AMF, IDAMAX, IZAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSTEBZ, F02EAZ, P01ABW, P01ABX, P01ABY,
     *                  ZDSCAL, ZHETRD, ZSCAL, ZSTEIN, ZSWAP, ZUNMTR
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCONJG, DIMAG, MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters
C
      WANTZ = (JOB.EQ.'V') .OR. (JOB.EQ.'v')
      UPPER = (UPLO.EQ.'U') .OR. (UPLO.EQ.'u')
      VALEIG = (RANGE.EQ.'V') .OR. (RANGE.EQ.'v')
C
      IERR = 0
      M = 0
      NREC = 0
      IF ( .NOT. WANTZ .AND. JOB.NE.'N' .AND. JOB.NE.'n')
     *    CALL P01ABW(JOB,'JOB',IFAIL,IERR,SRNAME)
      IF ( .NOT. VALEIG .AND. RANGE.NE.'I' .AND. RANGE.NE.'i')
     *    CALL P01ABW(RANGE,'RANGE',IFAIL,IERR,SRNAME)
      IF ( .NOT. UPPER .AND. UPLO.NE.'L' .AND. UPLO.NE.'l')
     *    CALL P01ABW(UPLO,'UPLO',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.MAX(1,N)) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (VALEIG) THEN
         IF (WU.LE.WL) CALL P01ABX(WU,'WU',IFAIL,IERR,SRNAME)
      ELSE
         IF (IL.LT.1 .AND. N.GT.0) CALL P01ABY(IL,'IL',IFAIL,IERR,
     *                                  SRNAME)
         IF ((IU.GT.N .OR. IU.LT.IL) .AND. N.GT.0) CALL P01ABY(IU,'IU',
     *       IFAIL,IERR,SRNAME)
      END IF
      IF (MEST.LT.1) CALL P01ABY(MEST,'MEST',IFAIL,IERR,SRNAME)
      IF (LDZ.LT.1 .OR. (WANTZ .AND. LDZ.LT.N)) CALL P01ABY(LDZ,'LDZ',
     *    IFAIL,IERR,SRNAME)
      IF (LWORK.LT.MAX(1,2*N)) CALL P01ABY(LWORK,'LWORK',IFAIL,IERR,
     *                              SRNAME)
C
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IERR = 1
         NREC = 1
         GO TO 260
      END IF
C
C     Check that the diagonal elements of A are real
C
      DO 20 I = 1, N
         IF (DIMAG(A(I,I)).NE.ZERO) THEN
            WRITE (REC,FMT=99995) I
            IERR = 5
            NREC = 1
            GO TO 260
         END IF
   20 CONTINUE
C
C     Get machine constants
C
      RMAX = ONE/SQRT(SQRT(X02AMF()))
C
C     Scale matrix so that maximum magnitude of an element lies
C     in the range ONE to RMAX (magnitude is defined by
C     abs(dble(z)) + abs(dimag(z))
C
      DO 40 J = 1, N
         IF (UPPER) THEN
            IMAX = IZAMAX(J,A(1,J),1)
         ELSE
            IMAX = J - 1 + IZAMAX(N-J+1,A(J,J),1)
         END IF
         RWORK(J) = ABS(DBLE(A(IMAX,J))) + ABS(DIMAG(A(IMAX,J)))
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
      ID = 1
      IE = ID + N
      IRWRK = IE + N
      ITAU = 1
      IWRK = ITAU + N
      LWRK = LWORK - IWRK + 1
C
C     Reduce symmetric matrix to tridiagonal form
C
      CALL ZHETRD(UPLO,N,A,LDA,RWORK(ID),RWORK(IE),WORK(ITAU),WORK(IWRK)
     *            ,LWRK,INFO)
C
      IF (VALEIG) THEN
         IF (SCALE) THEN
            WLT = WL*SIGMA
            WUT = WU*SIGMA
         ELSE
            WLT = WL
            WUT = WU
         END IF
      END IF
C
C     Compute eigenvalues by bisection
C
      IF (WANTZ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      IBL = 1
      ISP = IBL + N
      IWO = ISP + N
      ABSTOL = ZERO
      CALL DSTEBZ(RANGE,ORDER,N,WLT,WUT,IL,IU,ABSTOL,RWORK(ID),RWORK(IE)
     *            ,M,NSPLIT,W,IWORK(IBL),IWORK(ISP),RWORK(IRWRK),
     *            IWORK(IWO),IERR)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99998)
         IERR = 2
         NREC = 1
         GO TO 240
      END IF
C
      IF (WANTZ) THEN
         IF (M.GT.MEST) THEN
            WRITE (REC,FMT=99997) M, MEST
            IERR = 3
            NREC = 4
            GO TO 220
         END IF
C
C        Compute eigenvectors by inverse iteration.
C
         IFA = IWO + N
         CALL ZSTEIN(N,RWORK(ID),RWORK(IE),M,W,IWORK(IBL),IWORK(ISP),Z,
     *               LDZ,RWORK(IRWRK),IWORK(IWO),IWORK(IFA),IERR)
         IF (IERR.GT.0) THEN
            WRITE (REC,FMT=99996) IERR
            IERR = 4
            NREC = 3
C
C           Set eigenvectors which have not converged to zero
C
            DO 120 J = 1, IERR
               JJ = IWORK(IFA+J-1)
               DO 100 I = 1, N
                  Z(I,J) = ZERO
  100          CONTINUE
  120       CONTINUE
         END IF
C
C        Sort selected eigenvalues, and reorder eigenvectors
C        correspondingly.
C
         DO 160 J = 1, M - 1
            I = 0
            TMP1 = W(J)
            DO 140 JJ = J + 1, M
               IF (W(JJ).LT.TMP1) THEN
                  I = JJ
                  TMP1 = W(JJ)
               END IF
  140       CONTINUE
            IF (I.NE.0) THEN
C
C              swap W(I) and W(J) and their corresponding eigenvectors
C
               W(I) = W(J)
               W(J) = TMP1
               CALL ZSWAP(N,Z(1,I),1,Z(1,J),1)
            END IF
  160    CONTINUE
C
C        Apply unitary matrix used in reduction to tridiagonal
C        form to eigenvectors returned by ZSTEIN.
C
         CALL ZUNMTR('Left',UPLO,'No transpose',N,M,A,LDA,WORK(ITAU),Z,
     *               LDZ,WORK(IWRK),LWRK,INFO)
C
C        Normalize eigenvectors so that element of largest absolute
C        value is positive.
C
         DO 200 J = 1, M
            DO 180 I = 1, N
               RWORK(I) = ABS(Z(I,J))
  180       CONTINUE
            JJ = IDAMAX(N,RWORK,1)
            TEMP = Z(JJ,J)/RWORK(JJ)
            CALL ZSCAL(N,DCONJG(TEMP),Z(1,J),1)
            Z(JJ,J) = RWORK(JJ)
  200    CONTINUE
      END IF
C
  220 CONTINUE
      IF (SCALE) THEN
C
C        Rescale eigenvalues
C
 
         CALL DSCAL(M,ONE/SIGMA,W,1)
      END IF
C
  240 CONTINUE
      IF (SCALE) THEN
C
C        Rescale tridiagonal matrix in A
C
         CALL ZDSCAL(N,ONE/SIGMA,A,LDA+1)
         IF (N.GT.1) THEN
            IF (UPPER) THEN
               CALL ZDSCAL(N-1,ONE/SIGMA,A(1,2),LDA+1)
            ELSE
               CALL ZDSCAL(N-1,ONE/SIGMA,A(2,1),LDA+1)
            END IF
         END IF
      END IF
C
  260 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** The input parameters contained ',I2,' error(s)')
99998 FORMAT (' ** The bisection algorithm failed to compute all the s',
     *       'pecified eigenvalues')
99997 FORMAT (' ** There are more than MEST eigenvalues in the specifi',
     *       'ed range.',/' ** M (number of eigenvalues in range) =',I6,
     *       '  MEST =',I6,/' ** No eigenvectors have been computed.',
     *       /' ** Rerun with 2nd dimension of Z = MEST >= M.')
99996 FORMAT (' ** Inverse iteration failed to compute all the specifi',
     *       'ed eigenvectors.',/' ** The number of eigenvectors which',
     *       ' failed to converge is',I6,/' ** The corresponding colum',
     *       'ns of Z are set to zero.')
99995 FORMAT (' ** A(i,i) has a non-zero imaginary part when i =',I6)
      END
