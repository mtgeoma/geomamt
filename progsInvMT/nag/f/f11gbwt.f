      SUBROUTINE F11GBW(SIGCMX,ITS,TALPHA,TBETA,D,E2,SIGMAX,SIGERR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBW - Compute the largest and smallest eigenvalues of T
C               (Symmetric iterative solver suite)
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGERR, SIGMAX, TALPHA, TBETA
      INTEGER           ITS, SIGCMX
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E2(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, EPSLON, MU1, MU1L, MUN, MUNH, SIGMX0,
     *                  SIGMX1, SIGTOL, TA1, TA2, TB1, TB2, TNORM, TOL,
     *                  X
      INTEGER           I, K, K1, K2, MAXITS, N, NI
      LOGICAL           SIGCMP
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           F11GBT
      EXTERNAL          X02AJF, F11GBT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, LOG, MAX, MIN
C     .. Save statement ..
      SAVE              EPS, EPSLON, MU1, MU1L, MUN, MUNH, SIGCMP,
     *                  SIGMX0, SIGMX1, SIGTOL, TA2, TB2, TNORM, MAXITS
C     .. Executable Statements ..
C
C     First call to F11GBW
C
      IF (SIGCMX.LT.0) THEN
         MAXITS = ITS
         ITS = 0
         SIGCMP = SIGCMX .LE. -2
         SIGCMX = -SIGCMX
         TA2 = TALPHA
         TB2 = ZERO
         IF (SIGCMP) THEN
            SIGTOL = TBETA
            TNORM = ABS(TA2)
            SIGMAX = TNORM
            SIGMX0 = -SIGMAX
            SIGMX1 = -SIGMAX
            SIGERR = ONE
            EPSLON = X02AJF()
            EPS = MAX(EPSLON,SIGTOL/(TWO*TWO))
CC            PRINT *, 'EPS ... ', EPS
            EPS = MAX(EPSLON,SIGTOL/(TWO**4))
CC            PRINT *, 'Enter EPS -- old EPS ... ', EPS
CC            READ *, EPS
            D(1) = TALPHA
            MU1L = TALPHA
            MU1 = TALPHA
            MUN = TALPHA
            MUNH = TALPHA
            TOL = EPS*TNORM
         ELSE
            SIGMAX = ABS(TA2)
            SIGERR = ZERO
         END IF
C
C     Subsequent call to F11GBW
C
      ELSE
C
C        Compute sigma[1] by bisection
C
         IF (SIGCMP) THEN
            ITS = ITS + 1
            N = ITS + 1
            TA1 = TA2
            TA2 = TALPHA
            TB1 = TB2
            TB2 = ABS(TBETA)
            TNORM = MAX(TNORM,(ABS(TA1)+TB2+TB1),(ABS(TA2)+TB2))
            D(N) = TALPHA
            E2(ITS) = TBETA**2
            MU1L = MIN(MU1L,(TA1-TB1-TB2),(TA2-TB2))
            MUNH = MAX(MUNH,(TA1+TB1+TB2),(TA2+TB2))
            TOL = TNORM*EPS
C
C           Compute lambda[1], the smallest eigenvalue of T[k], where
C           MU1L <= lambda[1] <= MU1.
C           Skip when MU1L >= -sigma_1(T[k-1])
C
            TA1 = MU1L
            IF (MU1L.LE.-SIGMAX) THEN
               X = MU1
               DO 20 I = 1, N
                  K2 = F11GBT(N,EPSLON,D,E2,MU1)
                  IF (K2.GE.1) GO TO 40
                  TA1 = MU1
                  MU1 = X + (MUNH-X)*DBLE(I)/DBLE(N)
   20          CONTINUE
   40          CONTINUE
C
               NI = INT(LOG((MU1-MU1L)/TOL)/LOG(TWO)) + 1
C
               DO 60 I = 1, NI
                  X = (TA1+MU1)/TWO
                  K = F11GBT(N,EPSLON,D,E2,X)
                  IF (K.LE.0) THEN
                     TA1 = X
C
C                    Exit from the loop when
C                    lambda[1] >= TA1 >= -sigma_1(T[k-1])
C
                     IF (TA1.GE.-SIGMAX) GO TO 80
                  ELSE
                     MU1 = X
                     K2 = K
                  END IF
                  IF (((MU1-TA1).LE.TOL) .AND. (K2.LE.1)) GO TO 80
   60          CONTINUE
   80          CONTINUE
            END IF
C
C           Compute lambda[n], the largest eigenvalue of T[k], where
C           MUN <= lambda[n] <= MUNH.
C           Skip when MUNH <= sigma_1(T[k-1])
C
            TB1 = MUNH
            IF (MUNH.GE.SIGMAX) THEN
               X = MUN
               DO 100 I = 1, N
                  K1 = F11GBT(N,EPSLON,D,E2,MUN)
                  IF (K1.LE.ITS) GO TO 120
                  TB1 = MUN
                  MUN = X + (MU1L-X)*DBLE(I)/DBLE(N)
  100          CONTINUE
  120          CONTINUE
C
               NI = INT(LOG((MUNH-MUN)/TOL)/LOG(TWO)) + 1
CC               PRINT *, 'NI ... ', NI
C
               DO 140 I = 1, NI
                  X = (MUN+TB1)/TWO
                  K = F11GBT(N,EPSLON,D,E2,X)
                  IF (K.LE.ITS) THEN
                     K1 = K
                     MUN = X
                  ELSE
                     TB1 = X
C
C                    Exit from the loop when
C                    lambda[1] <= TB1 <= sigma_1(T[k-1])
C
CC                     IF (TB1.LE.SIGMAX) GO TO 80
                     IF (TB1.LT.SIGMAX) GO TO 160
                  END IF
                  IF (((TB1-MUN).LE.TOL) .AND. (K1.GE.ITS)) GO TO 160
  140          CONTINUE
  160          CONTINUE
CC               PRINT *, 'I (final) ... ', I
            END IF
C
C           Complete
C
            TA1 = MAX(ABS(MU1),ABS(MUN))
CC            PRINT *, 'ITS, TNORM, TA1 ... ', ITS, TNORM, TA1
CC            PRINT *, 'ITS, TA1 ... ', ITS, TA1
CC            PRINT *, 'SIGMX0, SIGMX1 ... ', SIGMX0, SIGMX1
            IF (TA1.GE.MAX(SIGMX0,SIGMX1)) THEN
               SIGMX1 = SIGMX0
               SIGMX0 = SIGMAX
               SIGMAX = TA1
               TA1 = MAX(ABS(TA1),ABS(TB1),SIGMAX)
               SIGERR = MAX(ABS(TA1-SIGMX0),ABS(TA1-SIGMX1))
               IF (SIGMAX.NE.ZERO) SIGERR = SIGERR/SIGMAX
               IF ((SIGERR.LE.SIGTOL) .OR. (ITS.GE.MAXITS)) SIGCMX = 0
            ELSE
               SIGCMX = -1
            END IF
C
C        Use ||T||_1 as an approximation
C
         ELSE
            TA1 = TA2
            TA2 = TALPHA
            TB1 = TB2
            TB2 = ABS(TBETA)
            SIGMAX = MAX(SIGMAX,(ABS(TA1)+TB2+TB1),(ABS(TA2)+TB2))
         END IF
C
      END IF
C
C     End of function F11GBW
C
      RETURN
      END
