      SUBROUTINE G03FCF(TYPE,N,NDIM,D,X,LDX,STRESS,DFIT,ITER,IOPT,WK,
     *                  IWK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Non-metric multidimensional scaling
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03FCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STRESS
      INTEGER           IFAIL, IOPT, ITER, LDX, N, NDIM
      CHARACTER         TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N*(N-1)/2), DFIT(2*N*(N-1)), WK(15*N*NDIM),
     *                  X(LDX,NDIM)
      INTEGER           IWK(N*(N-1)/2+N*NDIM+5)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS
      INTEGER           I, IERROR, IFAULT, IRANKO, J, NIT, NM, NMISS,
     *                  NN, NREC, NZERO
      CHARACTER*22      OPTION
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04DGF, E04DKF, F06EFF, G03FCZ, M01DAF, M01EAF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (NDIM.LT.1) THEN
         WRITE (P01REC,FMT=99999) NDIM
      ELSE IF (N.LE.NDIM) THEN
         WRITE (P01REC,FMT=99998) N, NDIM
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC,FMT=99997) LDX, N
      ELSE IF (TYPE.NE.'T' .AND. TYPE.NE.'t' .AND. TYPE.NE.'S' .AND.
     *         TYPE.NE.'s') THEN
         WRITE (P01REC,FMT=99996) TYPE
      ELSE
         IERROR = 2
         NN = N*(N-1)/2
         NMISS = 0
         NZERO = 0
         DO 20 I = 1, NN
            IF (D(I).LT.0.0D0) NMISS = NMISS + 1
            IF (D(I).EQ.0.0D0) NZERO = NZERO + 1
   20    CONTINUE
         IF (NMISS+NZERO.EQ.NN) THEN
            WRITE (P01REC,FMT=99995)
         ELSE
            IERROR = 0
            NM = N*NDIM
C
C     Put information into workspace
C
            IRANKO = NM + 5
            IWK(NM+2) = N
            IWK(NM+3) = NDIM
            IWK(NM+4) = NMISS
            IF (TYPE.EQ.'T' .OR. TYPE.EQ.'t') THEN
               IWK(NM+5) = 0
            ELSE IF (TYPE.EQ.'S' .OR. TYPE.EQ.'s') THEN
               IWK(NM+5) = 1
            END IF
C
C     Copy initial values into workspace
C
            DO 60 J = 1, NDIM
               DO 40 I = 1, N
                  WK((J-1)*N+I) = X(I,J)
   40          CONTINUE
   60       CONTINUE
C
C     Minimize squared stress function by conjugate gradient
C
C     Sort distances
C
            IFAULT = 0
            CALL M01DAF(D,1,NN,'A',IWK(IRANKO+1),IFAULT)
C
C     Set options
C
            IF (IOPT.GE.0) THEN
               CALL E04DKF('Nolist')
               CALL E04DKF('Defaults')
               CALL E04DKF('Nolist')
               IF (IOPT.EQ.0) THEN
                  CALL E04DKF('Optimality tolerance = 0.00001')
               ELSE
                  EPS = 10.0D0**(-IOPT)
                  WRITE (OPTION,FMT=99989) 'O=', EPS
                  CALL E04DKF(OPTION)
               END IF
               CALL E04DKF('Print level = -1')
               CALL E04DKF('Verify = No')
            END IF
            IF (ITER.EQ.0) THEN
               CALL E04DKF('Iters=50')
            ELSE IF (ITER.GT.0) THEN
               WRITE (OPTION,FMT=99990) 'Iters=', ITER
               CALL E04DKF(OPTION)
            ELSE IF (ITER.LT.0) THEN
               NIT = MAX(50,5*NM)
               WRITE (OPTION,FMT=99990) 'Iters=', NIT
               CALL E04DKF(OPTION)
            END IF
            NIT = 0
            IFAULT = 1
            CALL E04DGF(NM,G03FCZ,NIT,STRESS,WK(NM+1),WK,IWK,WK(2*NM+1),
     *                  IWK(NM+2),DFIT,IFAULT)
            IF (IFAULT.EQ.3) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99994)
            ELSE IF (IFAULT.EQ.6) THEN
               IERROR = 4
               WRITE (P01REC,FMT=99993)
            ELSE IF (IFAULT.EQ.8) THEN
               IERROR = 5
               WRITE (P01REC,FMT=99992)
            ELSE IF (IFAULT.NE.0) THEN
               IERROR = 6
               WRITE (P01REC,FMT=99991)
            END IF
            IWK(1) = NIT
C
C     Decompress matrix X
C
            DO 100 J = 1, NDIM
               DO 80 I = 1, N
                  X(I,J) = WK((J-1)*N+I)
   80          CONTINUE
  100       CONTINUE
C
C     Reset ordered distances
C
            IFAULT = 0
            CALL F06EFF(NN,DFIT,1,DFIT(NN+1),1)
            CALL M01EAF(DFIT(NN+1),1,NN,IWK(IRANKO+1),IFAULT)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, NDIM .lt. 1 : NDIM = ',I16)
99998 FORMAT (' ** On entry, N .le. NDIM: N = ',I16,' NDIM = ',I16)
99997 FORMAT (' ** On entry, LDX .lt. N: LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, TYPE is not valid: TYPE = ',A1)
99995 FORMAT (' ** On entry, all the elements of D <= 0.0')
99994 FORMAT (' ** Optimization has not converged in given number ',
     *       'of iterations')
99993 FORMAT (' ** No acceptable solution found')
99992 FORMAT (' ** Optimization is not possible from initial configura',
     *       'tion')
99991 FORMAT (' ** Optimization has failed')
99990 FORMAT (A,I16)
99989 FORMAT (A,E13.5)
      END
