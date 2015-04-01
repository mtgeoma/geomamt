      SUBROUTINE G07BBF(METHOD,N,X,XC,IC,XMU,XSIG,TOL,MAXIT,SEXMU,
     *                  SEXSIG,CORR,DEV,NOBS,NIT,WK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     PROGRAM COMPUTES MAXIMUM-LIKELIHOOD ESTIMATES FOR PARAMETERS OF
C     THE NORMAL DISTRIBUTION FROM GROUPED AND/OR CENSORED DATA USING
C     EITHER A NEWTON-RAPHSON OR EXPECTATION-MAXIMIZATION APPROACH.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CORR, DEV, SEXMU, SEXSIG, TOL, XMU, XSIG
      INTEGER           IFAIL, MAXIT, N, NIT
      CHARACTER         METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  WK(2*N), X(N), XC(N)
      INTEGER           IC(N), NOBS(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  L11, L12, L22, LL, P, PE, PI, PP, RR, SS, SUM,
     *                  SUM2, SUMG, T1, T2, TEMP, TOLM
      INTEGER           CCODE, I, IERROR, IFAIL2, MAXITS, N1, N2, N3,
     *                  NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X01AAF, X02AJF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          S15ABF, X01AAF, X02AJF, X02ALF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G07BBY, G07BBZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG, DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IFAIL2 = 0
      CCODE = 0
      NOBS(1) = 0
      NOBS(2) = 0
      NOBS(3) = 0
      NOBS(4) = 0
      IERROR = 0
      PP = 0.0D0
      SUM = 0.0D0
      SUMG = 0.0D0
      SUM2 = 0.0D0
      PI = X01AAF(0.0D0)
C
      IF (X02AJF().LE.TOL .AND. TOL.LE.1.0D0) THEN
         TOLM = TOL
      ELSE IF (TOL.EQ.0.0D0) THEN
         TOLM = 0.000005D0
      ELSE
         IERROR = 1
         WRITE (P01REC,FMT=99999) TOL
         GO TO 120
      END IF
      IF (METHOD.NE.'N' .AND. METHOD.NE.'n' .AND. METHOD.NE.'E' .AND.
     *    METHOD.NE.'e') THEN
         IERROR = 1
         WRITE (P01REC,FMT=99996) METHOD
         GO TO 120
      END IF
      IF (MAXIT.LE.0) THEN
         MAXITS = 25
      ELSE
         MAXITS = MAXIT
      END IF
      IF (N.LT.2) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99991) N
      ELSE
         DO 20 I = 1, N
            IF (IC(I).EQ.1) THEN
               NOBS(1) = NOBS(1) + 1
            ELSE IF (IC(I).EQ.2) THEN
               NOBS(2) = NOBS(2) + 1
            ELSE IF (IC(I).EQ.3) THEN
               IF (X(I).NE.XC(I)) THEN
                  NOBS(3) = NOBS(3) + 1
                  SUMG = SUMG + 0.5D0*(X(I)+XC(I))
               ELSE
                  CCODE = CCODE + 1
               END IF
            ELSE IF (IC(I).EQ.0) THEN
               NOBS(4) = NOBS(4) + 1
               SUM = SUM + X(I)
               SUM2 = SUM2 + X(I)**2
            ELSE
               IERROR = 1
               WRITE (P01REC,FMT=99998) I, IC(I)
               GO TO 120
            END IF
   20    CONTINUE
         IF ((N-CCODE).LT.2) THEN
            IERROR = 1
            WRITE (P01REC,FMT=99997)
         ELSE
            IF (XSIG.LE.0.0D0) THEN
               IF (NOBS(4).GE.2) THEN
                  SS = 0.0D0
                  XMU = SUM/DBLE(NOBS(4))
                  DO 40 I = 1, N
                     IF (IC(I).EQ.0) SS = SS + (X(I)-XMU)**2
   40             CONTINUE
                  XSIG = SQRT(SS/(DBLE(NOBS(4))-1))
               ELSE IF (NOBS(3).GE.1) THEN
                  SS = 0.0D0
                  XMU = SUMG/DBLE(NOBS(3))
                  DO 60 I = 1, N
                     IF (IC(I).EQ.3) THEN
                        IF (X(I).NE.XC(I)) SS = SS + (0.5D0*(X(I)-XC(I))
     *                      -XMU)**2
                     END IF
   60             CONTINUE
                  XSIG = SQRT(SS/DBLE(NOBS(3)))
               ELSE
                  XMU = 0.0D0
                  XSIG = 1.0D0
               END IF
            END IF
            N1 = 1
            N2 = NOBS(1) + 1
            N3 = N2 + NOBS(2)
            DO 80 I = 1, N
               IF (IC(I).EQ.1) THEN
                  WK(N1) = X(I)
                  N1 = N1 + 1
               ELSE IF (IC(I).EQ.2) THEN
                  WK(N2) = X(I)
                  N2 = N2 + 1
               ELSE IF (IC(I).EQ.3) THEN
                  IF (X(I).LT.XC(I)) THEN
                     WK(N3) = X(I)
                     WK(N3+1) = XC(I)
                     N3 = N3 + 2
                  ELSE IF (X(I).GT.XC(I)) THEN
                     WK(N3) = XC(I)
                     WK(N3+1) = X(I)
                     N3 = N3 + 2
                  END IF
               END IF
   80       CONTINUE
            IF (METHOD.EQ.'E' .OR. METHOD.EQ.'e') THEN
               CALL G07BBY(XMU,XSIG,NOBS,TOLM,MAXITS,SUM,SUM2,L11,L12,
     *                     L22,NIT,WK,IERROR)
C
            ELSE IF (METHOD.EQ.'N' .OR. METHOD.EQ.'n') THEN
               CALL G07BBZ(XMU,XSIG,NOBS,TOLM,MAXITS,SUM,SUM2,L11,L12,
     *                     L22,NIT,WK,IERROR)
            END IF
            IF (IERROR.EQ.2) WRITE (P01REC,FMT=99995) MAXITS
            IF (IERROR.NE.3) THEN
               DO 100 I = 1, N
                  IF (IC(I).EQ.1) THEN
                     RR = (X(I)-XMU)/XSIG
                     P = 1.0D0 - S15ABF(RR,IFAIL2)
                     IF (P.GT.0.0D0) PP = PP + LOG(P)
                  ELSE IF (IC(I).EQ.2) THEN
                     RR = (X(I)-XMU)/XSIG
                     P = S15ABF(RR,IFAIL2)
                     IF (P.GT.0.0D0) PP = PP + LOG(P)
                  ELSE IF (IC(I).EQ.3) THEN
                     RR = (X(I)-XMU)/XSIG
                     LL = (XC(I)-XMU)/XSIG
                     P = S15ABF(LL,IFAIL2) - S15ABF(RR,IFAIL2)
                     IF (P.GT.0.0D0) PP = PP + LOG(P)
                  END IF
  100          CONTINUE
               PE = DBLE(NOBS(4))
               DEV = -PE*LOG(XSIG*SQRT(2.0D0*PI)) -
     *               (SUM2-2.0D0*XMU*SUM+PE*XMU*XMU)/(2.0D0*XSIG*XSIG)
               DEV = DEV + PP
               IF (ABS(L22).LT.1.0D0) THEN
                  TEMP = (L11*L22-L12**2)
               ELSE IF (ABS(L11).GT.X02ALF()/ABS(L22)) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99992)
                  GO TO 120
               ELSE IF (ABS(L12).LT.1.0D0) THEN
                  TEMP = (L11*L22-L12**2)
               ELSE IF (ABS(L12).GT.X02ALF()/ABS(L12)) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99992)
                  GO TO 120
               ELSE
                  TEMP = (L11*L22-L12**2)
               END IF
               IF (TEMP.EQ.0.0D0) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99992)
               ELSE
                  T1 = -L22/TEMP
                  T2 = -L11/TEMP
                  IF (T1.LE.0.0D0 .OR. T2.LE.0.0D0) THEN
                     IERROR = 4
                     WRITE (P01REC,FMT=99992)
                  ELSE
                     SEXMU = SQRT(T1)
                     SEXSIG = SQRT(T2)
                     CORR = (L12/TEMP)/(SEXMU*SEXSIG)
                  END IF
               END IF
            ELSE IF (METHOD.EQ.'N' .OR. METHOD.EQ.'n') THEN
               WRITE (P01REC,FMT=99993)
            ELSE
               WRITE (P01REC,FMT=99994)
            END IF
         END IF
      END IF
C
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, TOL is invalid : TOL = ',D13.5)
99998 FORMAT (' ** On entry, IC(',I16,') is invalid, it contains',I16)
99997 FORMAT (' ** On entry, effective number of observations .lt. 2 ')
99996 FORMAT (' ** On entry, METHOD is not a valid character : METHOD ',
     *       '= ',A1)
99995 FORMAT (' ** Method has not converged in ',I16,' iterations')
99994 FORMAT (' ** The EM process has failed')
99993 FORMAT (' ** Process has diverged')
99992 FORMAT (' ** Standard errors cannot be computed')
99991 FORMAT (' ** On entry, N.lt.2 : N = ',I16)
      END
