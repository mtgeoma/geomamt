      SUBROUTINE G12BAF(OFFSET,N,M,NS,Z,LDZ,ISZ,IP,T,IC,OMEGA,ISI,DEV,B,
     *                  SE,SC,COV,RES,ND,TP,SUR,NDMAX,TOL,MAXIT,IPRINT,
     *                  WK,IWK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C
C     PURPOSE:
C
C     This subroutine computes:
C         1. Maximum Likelihood Estimates for the parameters of
C             Cox's proportional harzard model using Marginal
C             Likelihood Approach.
C         2. Asymptotic standard deviate  and its P-value for 2-sided
C            test
C         3. Variance-Covariance Matrix of the ML estimates
C         4. Hazard function for each distinct failure time
C
C

C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G12BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DEV, TOL
      INTEGER           IFAIL, IP, IPRINT, LDZ, M, MAXIT, N, ND, NDMAX,
     *                  NS
      CHARACTER         OFFSET
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV(IP*(IP+1)/2), OMEGA(*), RES(N),
     *                  SC(IP), SE(IP), SUR(NDMAX,*), T(N), TP(NDMAX),
     *                  WK(IP*(IP+9)/2+N), Z(LDZ,M)
      INTEGER           IC(N), ISI(*), ISZ(M), IWK(2*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  TIME, TOLA
      INTEGER           I, IERROR, IFAULT, ISTRAT, ISWAP, ITER, J, K,
     *                  L1, L2, L3, L4, L5, LASTC, N0, NREC
      LOGICAL           NOCEN, NODEAD, OFFL
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G12BAZ, M01DAF, M01ZAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      TOLA = 1.0D0*X02AJF()
      IF (M.LT.1) THEN
         WRITE (P01REC,FMT=99999) M
      ELSE IF (N.LT.2) THEN
         WRITE (P01REC,FMT=99998) N
      ELSE IF (NS.LT.0) THEN
         WRITE (P01REC,FMT=99997) NS
      ELSE IF (LDZ.LT.N) THEN
         WRITE (P01REC,FMT=99996) LDZ
      ELSE IF (MAXIT.LT.0) THEN
         WRITE (P01REC,FMT=99993) MAXIT
      ELSE IF (TOL.LT.TOLA) THEN
         WRITE (P01REC,FMT=99991) TOL
      ELSE IF (OFFSET.NE.'Y' .AND. OFFSET.NE.'y' .AND. OFFSET.NE.
     *         'N' .AND. OFFSET.NE.'n') THEN
         WRITE (P01REC,FMT=99992) OFFSET
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.1) GO TO 260
C
C     Check ISZ and IP
C
      J = 0
      DO 20 I = 1, M
         IF (ISZ(I).LT.0) THEN
            IERROR = 2
            WRITE (P01REC,FMT=99985) I
            GO TO 260
         ELSE IF (ISZ(I).GT.0) THEN
            J = J + 1
         END IF
   20 CONTINUE
      IF (J.NE.IP) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99984)
         GO TO 260
      END IF
C
C     Check strata indicator
C
      IF (NS.GT.0) THEN
         N0 = 0
         DO 40 I = 1, N
            J = ISI(I)
            IF (J.LT.0 .OR. J.GT.NS) THEN
               IERROR = 2
               NREC = 2
               WRITE (P01REC,FMT=99989) I, J, NS
               GO TO 260
            ELSE IF (J.GT.0) THEN
               N0 = N0 + 1
            END IF
   40    CONTINUE
      ELSE
         N0 = N
      END IF
      IF (N0.LE.IP) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99986)
         GO TO 260
      END IF
C
C     Index T and IC and find distict failure times
C
      IFAULT = 0
      CALL M01DAF(T,1,N,'A',IWK,IFAULT)
      CALL M01ZAF(IWK,1,N,IFAULT)
C
C     Compute distinct failure
C
      IF (NS.GT.0) THEN
         DO 60 J = 1, N
            I = IWK(J)
            IF (ISI(I).GT.0) GO TO 80
   60    CONTINUE
   80    CONTINUE
      ELSE
         ISTRAT = 1
         J = 1
         I = IWK(1)
      END IF
      TIME = T(I)
      IF (IC(I).EQ.0) THEN
         ND = 1
         NOCEN = .TRUE.
         NODEAD = .FALSE.
         IWK(N+1) = J
      ELSE IF (IC(I).EQ.1) THEN
         ND = 0
         NOCEN = .FALSE.
         NODEAD = .TRUE.
         LASTC = J
      ELSE
         IERROR = 2
         I = J
         WRITE (P01REC,FMT=99990) J, IC(J)
         GO TO 260
      END IF
      DO 100 I = J + 1, N
         K = IWK(I)
         IF (NS.GT.0) THEN
            ISTRAT = ISI(K)
         END IF
         IF (ISTRAT.GT.0) THEN
            IF (T(K).EQ.TIME .AND. IC(K).EQ.0 .AND. .NOT. NOCEN) THEN
               IF (NODEAD) THEN
                  ND = ND + 1
                  IWK(N+ND) = LASTC
                  NODEAD = .FALSE.
               END IF
               ISWAP = IWK(LASTC)
               IWK(LASTC) = K
               IWK(I) = ISWAP
               LASTC = LASTC + 1
            ELSE IF (T(K).EQ.TIME .AND. IC(K).EQ.1 .AND. NOCEN) THEN
               NOCEN = .FALSE.
               LASTC = I
            ELSE IF (T(K).GT.TIME .AND. IC(K).EQ.0) THEN
               ND = ND + 1
               IWK(N+ND) = I
               NOCEN = .TRUE.
               NODEAD = .FALSE.
            ELSE IF (T(K).GT.TIME .AND. IC(K).EQ.1) THEN
               NOCEN = .FALSE.
               NODEAD = .TRUE.
               LASTC = I
            ELSE IF (IC(K).NE.1 .AND. IC(K).NE.0) THEN
               IERROR = 2
               WRITE (P01REC,FMT=99990) K, IC(K)
               GO TO 260
            END IF
            TIME = T(K)
         END IF
  100 CONTINUE
      IF (ND.LT.1) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99995)
         GO TO 260
      ELSE IF (NDMAX.LT.ND) THEN
         IERROR = 2
         NREC = 2
         WRITE (P01REC,FMT=99994) NDMAX, ND
         GO TO 260
      END IF
      DO 120 I = 1, ND
         K = IWK(N+I)
         TP(I) = T(IWK(K))
  120 CONTINUE
C
C    Compute means
C
      DO 140 J = 1, IP
         WK(J) = 0.0D0
  140 CONTINUE
      IF (NS.EQ.0) THEN
         DO 180 I = 1, N
            K = 0
            DO 160 J = 1, M
               IF (ISZ(J).NE.0) THEN
                  K = K + 1
                  WK(K) = WK(K) + Z(I,J)
               END IF
  160       CONTINUE
  180    CONTINUE
      ELSE
         DO 220 I = 1, N
            IF (ISI(I).NE.0) THEN
               K = 0
               DO 200 J = 1, M
                  IF (ISZ(J).NE.0) THEN
                     K = K + 1
                     WK(K) = WK(K) + Z(I,J)
                  END IF
  200          CONTINUE
            END IF
  220    CONTINUE
      END IF
      DO 240 J = 1, IP
         WK(J) = WK(J)/DBLE(N0)
  240 CONTINUE
C
C    Set up partitions of workspace
C
      L1 = IP + 1
      L2 = IP + L1
      L3 = IP + L2
      L4 = IP + L3
      L5 = IP*(IP+1)/2 + L4
      IF (OFFSET.EQ.'Y' .OR. OFFSET.EQ.'y') THEN
         OFFL = .TRUE.
      ELSE
         OFFL = .FALSE.
      END IF
C
      CALL G12BAZ(M,N,NS,Z,LDZ,WK,T,ISZ,IP,ISI,IC,OFFL,OMEGA,DEV,B,SE,
     *            SC,COV,SUR,NDMAX,IWK,IWK(N+1),WK(L1),WK(L2),WK(L3),
     *            WK(L4),WK(L5),RES,ND,MAXIT,IPRINT,TOL,ITER,IERROR)
C
      IWK(1) = ITER
      IF (IERROR.EQ.3) THEN
         WRITE (P01REC,FMT=99988)
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (P01REC,FMT=99982)
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (P01REC,FMT=99987) MAXIT
      ELSE IF (IERROR.EQ.6) THEN
         WRITE (P01REC,FMT=99983)
      END IF
  260 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, M .lt. 1: M = ',I16)
99998 FORMAT (' ** On entry, N .lt. 2: N = ',I16)
99997 FORMAT (' ** On entry, NS .lt. 0: NS = ',I16)
99996 FORMAT (' ** On entry, LDZ .lt. N: LDZ = ',I16)
99995 FORMAT (' ** All observations are censored.')
99994 FORMAT (' ** On entry, NDMAX too small, NDMAX = ',I16,/'        ',
     *       '                minimum value = ',I16)
99993 FORMAT (' ** On entry, MAXIT .lt. 0: MAXIT = ',I6)
99992 FORMAT (' ** On entry, OFFSET not valid: OFFSET = ',A1)
99991 FORMAT (' ** On entry, TOL .lt. 10*machine precision: TOL = ',
     *       E13.5)
99990 FORMAT (' ** On entry, IC must be 0 or 1: IC(',I16,') = ',I16)
99989 FORMAT (' ** On entry, ISI(',I16,') = ',I16,'.lt. 0',/'         ',
     *       '   or .gt. NS: NS = ',I16)
99988 FORMAT (' ** The matrix of second partial derivative is singular',
     *       '; computation abandoned.')
99987 FORMAT (' ** Convergence not achieved in ',I16,' iterations.')
99986 FORMAT (' ** On entry too few observations included in model.')
99985 FORMAT (' ** On entry the ',I16,' th value of ISZ .lt. 0.')
99984 FORMAT (' ** On entry there are not IP values of ISZ .gt. 0.')
99983 FORMAT (' ** Too many step halvings required. Convergence assume',
     *       'd.')
99982 FORMAT (' ** Overflow in calculations')
      END
