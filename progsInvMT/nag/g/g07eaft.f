      SUBROUTINE G07EAF(METHOD,N,X,CLEVEL,THETA,THETAL,THETAU,ESTCL,
     *                  WLOWER,WUPPER,WRK,IWRK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07EAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CLEVEL, ESTCL, THETA, THETAL, THETAU, WLOWER,
     *                  WUPPER
      INTEGER           IFAIL, N
      CHARACTER         METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(4*N), X(N)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, DIFF, EPS, ESTCL1, ESTCL2, GAMMA, SLOPE,
     *                  TOL, X1, X2
      INTEGER           I, IERROR, IF2, IND, IND1, IND2, INDH1, INDH2,
     *                  IX1, IX2, K1, K2, M, NREC
      LOGICAL           CHECKP, LOWER
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G07EAW, G07EAY, G07EAZ, M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MOD
C     .. Executable Statements ..
      IF (METHOD.NE.'E' .AND. METHOD.NE.'e' .AND. METHOD.NE.'A' .AND.
     *    METHOD.NE.'a') THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) METHOD
      ELSE IF (N.LT.2) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99998) N
      ELSE IF (CLEVEL.LE.0.0D0 .OR. CLEVEL.GE.1.0D0) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99997) CLEVEL
      ELSE
         NREC = 0
         IERROR = 0
         M = N*(N+1)/2
         DO 20 I = 1, N
            WRK(I) = X(I)
            IWRK(I) = I
   20    CONTINUE
         IF2 = 1
         CALL M01CAF(WRK,1,N,'A',IF2)
         IF (WRK(1).EQ.WRK(N)) THEN
C
C           All values identical - IFAIL = 2
C
            THETA = WRK(1)
            THETAL = THETA
            THETAU = THETA
            WLOWER = DBLE(M)
            WUPPER = 0.0D0
            IERROR = 2
            NREC = 3
            WRITE (P01REC,FMT=99996)
         ELSE
C
C           For the lower confidence limit we are in the upper tail of
C           the Wilcoxon distribution (non-increasing function).
C           Find INDH1 and INDH2 as the values of the Wilcoxon signed-
C           rank statistic that satisfy the required probabilities.
C
            LOWER = .TRUE.
            GAMMA = 0.5D0 - CLEVEL/2.0D0
            CALL G07EAZ(LOWER,GAMMA,N,INDH1,IWRK)
            WLOWER = DBLE(INDH1)
            IX1 = N - 2*INDH1/(N+1)
            IF (IX1.LT.1) IX1 = 1
            X1 = WRK(IX1)
C
            LOWER = .FALSE.
            GAMMA = 0.5D0 - CLEVEL/2.0D0
            CALL G07EAZ(LOWER,GAMMA,N,INDH2,IWRK)
            WUPPER = DBLE(INDH2)
            IX2 = N - 2*INDH2/(N+1)
            X2 = WRK(IX2)
            DIFF = X2 - X1
            SLOPE = (WUPPER-WLOWER)/DIFF
C            SLOPE1 = -INDH2/(WRK(N)-X2)
C
C           Use X1 and X2 as initial estimates of the limits to set
C           up the tolerance - 5 significant places relative to the
C           interval width.
C
            TOL = 0.25D-05*(X2-X1)
            EPS = 10.0D0*X02AJF()
            IF (TOL.LT.EPS) TOL = EPS
C
C           Find Lower Confidence Point
C
            CHECKP = .TRUE.
            LOWER = .TRUE.
            IF (METHOD.EQ.'E' .OR. METHOD.EQ.'e') THEN
               K1 = M - INDH1 + 1
               CALL G07EAW(WRK,N,K1,K1,LOWER,CHECKP,THETAL,ESTCL1,EPS,
     *                     IWRK,WRK(N+1))
            ELSE
               CALL G07EAY(WRK,N,INDH1,LOWER,CHECKP,SLOPE,THETAL,ESTCL1,
     *                     TOL,EPS,WRK(N+1),IERROR)
               IF (IERROR.EQ.3) THEN
                  NREC = 2
                  WRITE (P01REC,FMT=99995)
               END IF
            END IF
C
C           Find Upper Confidence Point
C
            LOWER = .FALSE.
            IF (METHOD.EQ.'E' .OR. METHOD.EQ.'e') THEN
               K2 = M - INDH2
               CALL G07EAW(WRK,N,K2,K2,LOWER,CHECKP,THETAU,ESTCL2,EPS,
     *                     IWRK,WRK(N+1))
            ELSE
               CALL G07EAY(WRK,N,INDH2,LOWER,CHECKP,SLOPE,THETAU,ESTCL2,
     *                     TOL,EPS,WRK(N+1),IERROR)
               IF (IERROR.EQ.3) THEN
                  NREC = 2
                  WRITE (P01REC,FMT=99994)
               END IF
            END IF
            ESTCL = 1.0D0 - (ESTCL1+ESTCL2)
C
C           Find estimate of theta.
C           First reset TOL and SLOPE using the final estimates of the
C           confidence interval endpoints.
C
            TOL = 0.25D-05*(THETAU-THETAL)
            IF (TOL.LT.EPS) TOL = EPS
            IF (THETAU.GT.THETAL) SLOPE = (WUPPER-WLOWER)
     *          /(THETAU-THETAL)
C
            CHECKP = .FALSE.
            IF (MOD(M,2).EQ.0) THEN
               IND1 = M/2
               IND2 = IND1 + 1
               IF (METHOD.EQ.'E' .OR. METHOD.EQ.'e') THEN
                  LOWER = .TRUE.
                  CALL G07EAW(WRK,N,IND1,IND2,LOWER,CHECKP,THETA,ESTCL1,
     *                        EPS,IWRK,WRK(N+1))
               ELSE
                  LOWER = .TRUE.
                  CALL G07EAY(WRK,N,IND1,LOWER,CHECKP,SLOPE,A1,ESTCL1,
     *                        TOL,EPS,WRK(N+1),IERROR)
                  IF (IERROR.EQ.3) THEN
                     NREC = 2
                     WRITE (P01REC,FMT=99993)
                  END IF
                  LOWER = .FALSE.
                  CALL G07EAY(WRK,N,IND1,LOWER,CHECKP,SLOPE,A2,ESTCL2,
     *                        TOL,EPS,WRK(N+1),IERROR)
                  IF (IERROR.EQ.3) THEN
                     NREC = 2
                     WRITE (P01REC,FMT=99993)
                  END IF
                  THETA = (A1+A2)/2.0D0
               END IF
            ELSE
               IND = (M+1)/2
               IF (METHOD.EQ.'E' .OR. METHOD.EQ.'e') THEN
                  LOWER = .TRUE.
                  CALL G07EAW(WRK,N,IND,IND,LOWER,CHECKP,THETA,ESTCL1,
     *                        EPS,IWRK,WRK(N+1))
               ELSE
                  LOWER = .TRUE.
                  CALL G07EAY(WRK,N,IND,LOWER,CHECKP,SLOPE,THETA,ESTCL1,
     *                        TOL,EPS,WRK(N+1),IERROR)
                  IF (IERROR.EQ.3) THEN
                     NREC = 2
                     WRITE (P01REC,FMT=99993)
                  END IF
               END IF
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
99999 FORMAT (1X,'** On entry, METHOD is an invalid character : METHOD',
     *       ' = ',A1)
99998 FORMAT (1X,'** On entry, N.lt.2 : N = ',I16)
99997 FORMAT (1X,'** On entry CLEVEL is out of range : CLEVEL = ',D13.5)
99996 FORMAT (1X,'** Not enough information to compute an interval est',
     *       'imate since the',/4X,'whole sample is identical. The com',
     *       'mon value is returned in THETA, THETAL',/4X,'and THETAU.')
99995 FORMAT (1X,'** Warning. The iterative procedure to find an estim',
     *       'ate of the lower ',/4X,'confidence point had not converg',
     *       'ed in 100 iterations.')
99994 FORMAT (1X,'** Warning. The iterative procedure to find an estim',
     *       'ate of the upper ',/4X,'confidence point had not converg',
     *       'ed in 100 iterations.')
99993 FORMAT (1X,'** Warning. The iterative procedure to find an estim',
     *       'ate of Theta ',/4X,'had not converged in 100 iterations.')
      END
