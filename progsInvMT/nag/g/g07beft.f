      SUBROUTINE G07BEF(CENS,N,X,IC,BETA,GAMMA,TOL,MAXIT,SEBETA,SEGAM,
     *                  CORR,DEV,NIT,WK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     PROGRAM FOR ESTIMATING PARAMETERS OF THE WEIBULL DISTRIBUTION.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07BEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BETA, CORR, DEV, GAMMA, SEBETA, SEGAM, TOL
      INTEGER           IFAIL, MAXIT, N, NIT
      CHARACTER         CENS
C     .. Array Arguments ..
      DOUBLE PRECISION  WK(N), X(N)
      INTEGER           IC(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, DG, DL, DXBETA, DXGAM, ETOL, L1, L11, L12,
     *                  L2, L22, RELBET, RELGAM, RN, RNE, S, SAFE, SLX,
     *                  SXG, SXGLX, SXGLX2, SXOM, T, TEMP, XMAX, ZMAX,
     *                  ZMIN
      INTEGER           I, IERROR, J, MAXITS, NE, NREC, QG, QL
      LOGICAL           CEN
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G07BEZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
      SAFE = -LOG(X02AMF())
C
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE
         IF (CENS.EQ.'C' .OR. CENS.EQ.'c') THEN
            CEN = .TRUE.
         ELSE IF (CENS.EQ.'N' .OR. CENS.EQ.'n') THEN
            CEN = .FALSE.
         ELSE
            IERROR = 1
            WRITE (P01REC,FMT=99998) CENS
            GO TO 260
         END IF
         IF (X02AJF().LE.TOL .AND. TOL.LE.1.0D0) THEN
            ETOL = TOL
         ELSE IF (TOL.EQ.0.0D0) THEN
            ETOL = 0.000005D0
         ELSE
            IERROR = 1
            WRITE (P01REC,FMT=99997) TOL
            GO TO 260
         END IF
         IF (MAXIT.LT.1) THEN
            MAXITS = 25
         ELSE
            MAXITS = MAXIT
         END IF
         NE = 0
         ZMAX = 0.0D0
         ZMIN = X(1)
         IF ( .NOT. CEN) NE = N
         DO 20 I = 1, N
            IF (X(I).LE.0.0D0) THEN
               GO TO 240
            ELSE
               IF (X(I).LT.ZMIN) ZMIN = X(I)
               IF (X(I).GT.ZMAX) ZMAX = X(I)
               IF (CEN) THEN
                  IF (IC(I).EQ.0) THEN
                     NE = NE + 1
                  ELSE IF (IC(I).NE.1) THEN
                     GO TO 220
                  END IF
               END IF
            END IF
   20    CONTINUE
         XMAX = ZMAX
         ZMAX = MAX(-LOG(ZMIN),LOG(ZMAX))
         SAFE = SAFE - 2.0D0*LOG(ZMAX)
         IF (NE.EQ.0) THEN
            IERROR = 3
            WRITE (P01REC,FMT=99995)
         ELSE
            IF (GAMMA.LE.0.0D0) THEN
               DO 40 I = 1, N
                  IF ( .NOT. CEN) THEN
                     WK(I) = X(I)
                  ELSE IF (IC(I).EQ.0) THEN
                     WK(I) = X(I)
                  ELSE
                     WK(I) = -X(I)
                  END IF
   40          CONTINUE
               CALL G07BEZ(N,NE,WK,GAMMA,IERROR)
               IF (IERROR.EQ.3) THEN
                  WRITE (P01REC,FMT=99989)
                  GO TO 260
               END IF
            END IF
C
C           SCALE DATA TO AVOID OVERFLOW
C
            SLX = 0.0D0
            SXOM = 0.0D0
            DO 60 I = 1, N
               AA = X(I)/XMAX
               SXOM = SXOM + (AA)**GAMMA
               WK(I) = LOG(X(I))
               IF ( .NOT. CEN) THEN
                  SLX = SLX + WK(I)
               ELSE IF (IC(I).EQ.0) THEN
                  SLX = SLX + WK(I)
               END IF
   60       CONTINUE
            RN = DBLE(N)
            RNE = DBLE(NE)
            BETA = LOG(RN) - LOG(SXOM) - GAMMA*LOG(XMAX)
C
C           START ITERATIONS
C
            QL = -1
            QG = -1
            NIT = 0
            DL = 0.0D0
            DG = 0.0D0
            DO 100 I = 1, MAXITS
               SXG = 0.0D0
               SXGLX = 0.0D0
               SXGLX2 = 0.0D0
               NIT = NIT + 1
               DO 80 J = 1, N
                  T = BETA + GAMMA*WK(J)
                  IF (ABS(T).GE.SAFE) THEN
                     GO TO 200
                  ELSE
                     S = EXP(T)
                     SXG = SXG + S
                     SXGLX = SXGLX + S*WK(J)
                     SXGLX2 = SXGLX2 + S*(WK(J)**2)
                  END IF
   80          CONTINUE
               L1 = RNE - SXG
               L2 = RNE/GAMMA + SLX - SXGLX
               L11 = -SXG
               L12 = -SXGLX
               L22 = -RNE/(GAMMA**2) - SXGLX2
               TEMP = (L11*L22-L12*L12)
               IF (TEMP.EQ.0.0D0) THEN
                  GO TO 180
               ELSE
                  DXBETA = (-L1*L22+L2*L12)/TEMP
                  DXGAM = (L1*L12-L2*L11)/TEMP
                  BETA = BETA + DXBETA
                  GAMMA = GAMMA + DXGAM
                  IF (GAMMA.LE.0.0D0) GAMMA = 0.5D0*(GAMMA-DXGAM)
                  RELBET = ABS(DXBETA/(1.0D0+ABS(BETA)))
                  RELGAM = ABS(DXGAM/(1.0D0+GAMMA))
                  IF (ABS(DXBETA).GT.DL) THEN
                     QL = QL + 1
                  ELSE
                     QL = 0
                  END IF
                  IF (ABS(DXGAM).GT.DG) THEN
                     QG = QG + 1
                  ELSE
                     QG = 0
                  END IF
                  DL = ABS(DXBETA)
                  DG = ABS(DXGAM)
                  IF (QL.GE.3 .OR. QG.GE.3) THEN
                     GO TO 160
                  ELSE IF (RELBET.LE.ETOL .AND. RELGAM.LE.ETOL) THEN
                     GO TO 120
                  END IF
               END IF
  100       CONTINUE
            IERROR = 4
            WRITE (P01REC,FMT=99990) MAXITS
            GO TO 260
  120       SEBETA = SQRT(-L22/TEMP)
            SEGAM = SQRT(-L11/TEMP)
            CORR = L12/(SQRT(L12*L22))
            DEV = 0.0D0
            DO 140 J = 1, N
               T = BETA + GAMMA*WK(J)
               IF (ABS(T).GE.SAFE) THEN
                  GO TO 200
               ELSE
                  DEV = DEV + EXP(T)
               END IF
  140       CONTINUE
            DEV = RNE*LOG(GAMMA) + RNE*BETA + (GAMMA-1.0D0)*SLX - DEV
            GO TO 260
  160       IERROR = 5
            WRITE (P01REC,FMT=99991)
            GO TO 260
  180       IERROR = 5
            WRITE (P01REC,FMT=99992)
            GO TO 260
  200       IERROR = 6
            WRITE (P01REC,FMT=99994)
         END IF
         GO TO 260
  220    IERROR = 2
         WRITE (P01REC,FMT=99993) I, IC(I)
         GO TO 260
  240    IERROR = 2
         WRITE (P01REC,FMT=99996) I, X(I)
      END IF
  260 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.1 : N =',I16)
99998 FORMAT (' ** On entry, CENS is not a valid character : CENS = ',
     *       A1)
99997 FORMAT (' ** On entry, TOL is invalid : TOL = ',D13.5)
99996 FORMAT (' ** On entry, the',I16,'th observation is .le. 0.0. X(I',
     *       ') =',D13.5)
99995 FORMAT (' ** On entry, there are no exact observations. ')
99994 FORMAT (' ** Potential overflow detected. ')
99993 FORMAT (' ** On entry, the',I16,' th IC was not valid. IC(I) = ',
     *       I16)
99992 FORMAT (' ** Hessian matrix is singular. ')
99991 FORMAT (' ** Iterations have diverged. ')
99990 FORMAT (' ** Iterations have failed to converge in ',I16,' itera',
     *       'tions.')
99989 FORMAT (' ** Unable to calculate initial values. ')
      END
