      SUBROUTINE G01JCF(A,MULT,RLAMDA,N,C,P,PDF,TOL,MAXIT,WRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01JCF')
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO, BETLIM
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,
     *                  BETLIM=1.8D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, P, PDF, TOL
      INTEGER           IFAIL, MAXIT, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), RLAMDA(N), WRK(N+2*MAXIT)
      INTEGER           MULT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A0, A0INV, ARG, BETA, DANS, HBETA, HOLD, HOLD2,
     *                  LANS, PANS, PREC, Q, Q1, SUM, SUM1, TOL2, UF,
     *                  UF1, UFLOW, UFLOW1, Z, ZL
      INTEGER           I, IERR, IFAIL2, II, K, M, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          S14BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MOD, DBLE, SQRT
C     .. Executable Statements ..
      NREC = 1
      P = ZERO
      PDF = ZERO
      IERR = 0
      IF (N.LT.1) THEN
         IERR = 1
         WRITE (REC,FMT=99999) N
      ELSE IF (MAXIT.LT.1) THEN
         IERR = 1
         WRITE (REC,FMT=99992) MAXIT
      ELSE IF (C.LT.ZERO) THEN
         IERR = 1
         WRITE (REC,FMT=99998) C
      ELSE IF (C.NE.ZERO) THEN
         PREC = X02AJF()*10.0D0
         IF (TOL.GT.PREC .AND. TOL.LT.1.0D0) PREC = TOL
         UFLOW1 = X02AMF()
         UFLOW = LOG(UFLOW1)
         SUM = A(1)
         BETA = SUM
         NREC = 2
         DO 20 I = 1, N
            HOLD = A(I)
            IF (HOLD.LE.ZERO) THEN
               IERR = 2
               WRITE (REC,FMT=99997) I, A(I)
               GO TO 160
            ELSE IF (MULT(I).LT.1) THEN
               IERR = 2
               WRITE (REC,FMT=99996) I, MULT(I)
               GO TO 160
            ELSE IF (RLAMDA(I).LT.ZERO) THEN
               IERR = 2
               WRITE (REC,FMT=99995) I, RLAMDA(I)
               GO TO 160
            ELSE IF (BETA.GT.HOLD) THEN
               BETA = HOLD
            ELSE IF (SUM.LT.HOLD) THEN
               SUM = HOLD
            END IF
   20    CONTINUE
         NREC = 1
         HBETA = BETA
         BETA = TWO/(ONE/BETA+ONE/SUM)
         IF (BETA.GT.BETLIM*HBETA) BETA = HBETA
         K = 0
         SUM = ZERO
         SUM1 = ZERO
         UF = ONE
         DO 40 I = 1, N
            HOLD = BETA/A(I)
            IF (ABS(ONE-HOLD).NE.ZERO .AND. ABS(ONE-HOLD).LT.UF)
     *          UF = ABS(ONE-HOLD)
            SUM = SUM + MULT(I)*LOG(HOLD)
            SUM1 = SUM1 + RLAMDA(I)
            K = K + MULT(I)
            WRK(I) = ONE
   40    CONTINUE
         UF1 = UFLOW1/UF
         ARG = HALF*(SUM-SUM1)
         IF (ARG.LT.UFLOW) THEN
            IERR = 5
            WRITE (REC,FMT=99993)
         ELSE
            A0 = EXP(ARG)
            Z = C/BETA
            ZL = LOG(Z)
            IF (MOD(K,2).EQ.0) THEN
               I = 2
               LANS = -HALF*Z
            ELSE
               I = 1
               LANS = -HALF*(Z+ZL) - LOG(SQRT(X01AAF(Q)/TWO))
            END IF
            K = K - 2
            DO 60 II = I, K, 2
               LANS = LANS + ZL - LOG(DBLE(II))
   60       CONTINUE
            IF (LANS.LT.UFLOW) THEN
               DANS = ZERO
            ELSE
               DANS = EXP(LANS)
            END IF
            IFAIL2 = 1
            CALL S14BAF((K+2)/TWO,Z/TWO,PREC,P,Q1,IFAIL2)
            IF (IFAIL2.NE.0) THEN
               IERR = 3
               WRITE (REC,FMT=99991)
               GO TO 160
            END IF
            PDF = DANS
            PANS = P
            TOL2 = PREC/A0
            A0INV = ONE/A0
            SUM = A0INV - ONE
            DO 120 M = 1, MAXIT
               SUM1 = ZERO
               DO 80 I = 1, N
                  HOLD = WRK(I)
                  IF (ABS(HOLD).GT.UF1) THEN
                     HOLD2 = HOLD - (HOLD*BETA)/A(I)
                  ELSE
                     HOLD2 = ZERO
                  END IF
                  WRK(I) = HOLD2
                  SUM1 = SUM1 + HOLD2*MULT(I) + M*RLAMDA(I)*(HOLD-HOLD2)
   80          CONTINUE
               SUM1 = HALF*SUM1
               WRK(N+MAXIT+M) = SUM1
               DO 100 I = M - 1, 1, -1
                  SUM1 = SUM1 + WRK(N+MAXIT+I)*WRK(N+M-I)
  100          CONTINUE
               SUM1 = SUM1/M
               WRK(N+M) = SUM1
               K = K + 2
               SUM = SUM - SUM1
               LANS = LANS + ZL - LOG(DBLE(K))
               IF (LANS.LT.UFLOW) THEN
                  DANS = ZERO
               ELSE IF (SUM1.EQ.ZERO) THEN
                  DANS = EXP(LANS)
                  PANS = PANS - DANS
               ELSE IF (LANS+LOG(ABS(SUM1)).LT.UFLOW) THEN
                  DANS = EXP(LANS)
                  PANS = PANS - DANS
               ELSE
                  DANS = EXP(LANS)
                  PANS = PANS - DANS
                  PDF = PDF + DANS*SUM1
               END IF
               SUM1 = PANS*SUM1
               P = P + SUM1
               IF (ABS(PANS*SUM).LT.TOL2 .AND. ABS(SUM1).LT.TOL2)
     *             GO TO 140
  120       CONTINUE
            IERR = 4
            WRITE (REC,FMT=99994) MAXIT
  140       P = A0*P
            PDF = A0*PDF/(BETA+BETA)
            IF (P.LT.ZERO .OR. P.GT.ONE) THEN
               IERR = 5
               WRITE (REC,FMT=99993)
               IF (P.LT.ZERO) THEN
                  P = ZERO
               ELSE
                  P = ONE
               END IF
               PDF = ZERO
            END IF
         END IF
      END IF
  160 IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
99999 FORMAT (1X,'** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (1X,'** On entry, C.lt.0.0 : C = ',D13.5)
99997 FORMAT (1X,'** On entry, A has an element .le. 0.0 ',/12X,': A(',
     *       I16,') = ',D13.5)
99996 FORMAT (1X,'** On entry, MULT has an element .lt. 1',/12X,': MUL',
     *       'T(',I16,') = ',I16)
99995 FORMAT (1X,'** On entry, RLAMDA has an element .lt. 0.0',/12X,
     *       ': RLAMDA(',I16,') = ',D13.5)
99994 FORMAT (1X,'** The required accuracy could not be met in ',I16,
     *       ' iterations.')
99993 FORMAT (1X,'** Calculated probability at boundary')
99992 FORMAT (1X,'** On entry, MAXIT.lt.1: MAXIT = ',I16)
99991 FORMAT (1X,'** The central Chi square has failed to converge')
      END
