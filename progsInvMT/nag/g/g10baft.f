      SUBROUTINE G10BAF(N,X,WINDOW,SLO,SHI,NS,SMOOTH,T,USEFFT,FFT,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C       Algorithm AS 176 APPL. STATIST. (1982) VOL. 31, NO.1
C
C       Find density estimate by kernel method using Gaussian
C       kernel. The interval on which the estimate is evaluated
C       has end points SLO and SHI. If USEFT is not zero
C       then it is assumed that the routine has been
C       called before with the same data and end points
C       and that the array FT has not been altered.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G10BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SHI, SLO, WINDOW
      INTEGER           IFAIL, N, NS
      LOGICAL           USEFFT
C     .. Array Arguments ..
      DOUBLE PRECISION  FFT(NS), SMOOTH(NS), T(NS), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AINC, BIG, DL, DM, FAC, FAC1, HW, L, M, ONE, RJ,
     *                  SHSAVE, SLO1, SLSAVE, STEP, TEMP, WINC, WT, WW,
     *                  ZERO
      INTEGER           I, IERROR, IFAIL2, J, J1, J2, J2LO, JHI, JJ,
     *                  JMAX, KK, NREC, NS2, NSAVE, NSSAVE
      LOGICAL           SAVED
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06EAF, C06EBF, C06GBF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, EXP, INT, LOG, MIN, SQRT
C     .. Save statement ..
      SAVE              SAVED, NSAVE, NSSAVE, SLSAVE, SHSAVE, L, M
C     .. Data statements ..
      DATA              ZERO, ONE/0.0D0, 1.0D0/
      DATA              SAVED/.FALSE./
C     .. Executable Statements ..
C
C       Initialize and check for valid parameter values
C
      IERROR = 0
      NREC = 1
      IF (N.LE.0) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) N
      ELSE IF (NS.LT.2) THEN
         IERROR = 1
         WRITE (REC,FMT=99998) NS
      ELSE IF (SHI.LE.SLO) THEN
         IERROR = 1
         WRITE (REC,FMT=99997) SHI, SLO
      ELSE IF (WINDOW.LE.0.0D0) THEN
         IERROR = 1
         WRITE (REC,FMT=99996) WINDOW
      ELSE
C
         BIG = -LOG(X02AMF())
         STEP = (SHI-SLO)/DBLE(NS)
         AINC = ONE/(DBLE(N)*STEP)
         NS2 = NS/2
         HW = WINDOW/STEP
         FAC1 = (X01AAF(TEMP)*HW/DBLE(NS))**2
         IF (USEFFT) THEN
C
C           Check scalar inputs
C
            IF ( .NOT. SAVED) THEN
               IERROR = 2
               NREC = 2
               WRITE (REC,FMT=99993)
               GO TO 120
            END IF
            IF ((N.NE.NSAVE) .OR. (NS.NE.NSSAVE) .OR. (SLO.NE.SLSAVE)
     *           .OR. (SHI.NE.SHSAVE)) THEN
               IERROR = 2
               NREC = 2
               WRITE (REC,FMT=99992)
               GO TO 120
            END IF
         ELSE
C
C           Discretize the data
C
            SLO1 = SLO - STEP*0.5D0
            DO 20 J = 1, NS
               FFT(J) = ZERO
   20       CONTINUE
            M = X(1)
            L = X(1)
            DO 40 I = 1, N
               TEMP = X(I)
               WT = (TEMP-SLO1)/STEP
               JJ = INT(WT)
               IF (JJ.LT.1 .OR. JJ.GT.NS) GO TO 40
               WT = WT - DBLE(JJ)
               WINC = WT*AINC
               KK = JJ + 1
               IF (JJ.EQ.NS) KK = 1
               FFT(JJ) = FFT(JJ) + AINC - WINC
               FFT(KK) = FFT(KK) + WINC
               IF (TEMP.LT.M) THEN
                  M = TEMP
               ELSE IF (TEMP.GT.L) THEN
                  L = TEMP
               END IF
   40       CONTINUE
C
C           Transform to find FFT
C
            IFAIL2 = 1
            CALL C06EAF(FFT,NS,IFAIL2)
            IF (IFAIL2.NE.0) THEN
               IERROR = 3
               WRITE (REC,FMT=99994)
               GO TO 120
            END IF
C
C           Save scalar inputs
C
            SAVED = .TRUE.
            NSAVE = N
            NSSAVE = NS
            SLSAVE = SLO
            SHSAVE = SHI
         END IF
C
C           Check sufficiency of bounds
C
         DM = M - 3.0D0*WINDOW
         DL = L + 3.0D0*WINDOW
         IF ((SLO.GT.DM) .OR. (SHI.LT.DL)) THEN
            IERROR = 4
            WRITE (REC,FMT=99995)
         END IF
C
C        Find transform of density estimate
C
         JHI = SQRT(BIG/FAC1)
         JMAX = MIN(NS2-1,JHI)
         SMOOTH(1) = FFT(1)
         RJ = ZERO
         DO 60 J = 1, JMAX
            RJ = RJ + ONE
            FAC = EXP(-FAC1*RJ*RJ)
            J1 = J + 1
            J2 = NS - J + 1
            SMOOTH(J1) = FAC*FFT(J1)
            SMOOTH(J2) = FAC*FFT(J2)
   60    CONTINUE
C
C        Cope with underflow by setting tail of transform to zero
C
         WW = JHI + 1 - NS2
         IF (WW.LT.0.0D0) THEN
            J2LO = JHI + 2
            DO 80 J1 = J2LO, NS2
               J2 = NS - J1 + 2
               SMOOTH(J1) = ZERO
               SMOOTH(J2) = ZERO
   80       CONTINUE
            SMOOTH(NS2+1) = ZERO
         ELSE IF (WW.EQ.0.0D0) THEN
            SMOOTH(NS2+1) = ZERO
         ELSE IF (WW.GT.0.0D0) THEN
            SMOOTH(NS2+1) = EXP(-FAC1*DBLE(NS2)**2)*FFT(NS2+1)
         END IF
C
C        Invert Fourier transform of SMOOTH to get estimate
C        and estimate negative density values
C
         CALL C06GBF(SMOOTH,NS,IFAIL2)
C
         CALL C06EBF(SMOOTH,NS,IFAIL2)
C
         DO 100 J = 1, NS
            IF (SMOOTH(J).LT.ZERO) SMOOTH(J) = ZERO
            T(J) = SLO + (J-0.5D0)*(SHI-SLO)/NS
  100    CONTINUE
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.le.0: N = ',I16)
99998 FORMAT (' ** On entry, NS.lt.2: NS = ',I16)
99997 FORMAT (' ** On entry, SHI.le.SLO: SHI = ',E13.5,' SLO = ',E13.5)
99996 FORMAT (' ** On entry, WINDOW.le.0.0: WINDOW = ',E13.5)
99995 FORMAT (' ** On entry, not enough space allowed at bounds.')
99994 FORMAT (' ** On entry, at least one prime factor of NS is greate',
     *       'r than 19.')
99993 FORMAT (' ** On entry USEFFT = .true. and there has not been a p',
     *       'revious call to',/'    G10BAF with USEFFT = .false.')
99992 FORMAT (' ** On entry USEFFT = .true. and some of the scalar inp',
     *       'uts have been altered',/'       since the previous call ',
     *       'to G10BAF with USEFFT = .false.')
      END
