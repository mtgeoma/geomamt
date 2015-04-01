      SUBROUTINE C05PCZ(FCN,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,DIAG,MODE,
     *                  FACTOR,NPRINT,INFO,NFEV,NJEV,R,LR,QTF,WA1,WA2,
     *                  WA3,WA4)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05PCZ (based on MINPACK routine HYBRJ )
C
C     The purpose of C05PCZ is to find a zero of a system of
C     N nonlinear functions in N variables by a modification
C     of the POWELL hybrid method. The user must provide a
C     subroutine which calculates the functions and the Jacobian.
C
C     **********
C
C     Revised to call BLAS.
C     P.J.D. Mayes, J.J. Du Croz, NAG Central Office, September 1987.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, P1, P5, P001, P0001, ZERO
      PARAMETER         (ONE=1.0D0,P1=0.1D0,P5=0.5D0,P001=0.001D0,
     *                  P0001=0.0001D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FACTOR, XTOL
      INTEGER           INFO, LDFJAC, LR, MAXFEV, MODE, N, NFEV, NJEV,
     *                  NPRINT
C     .. Array Arguments ..
      DOUBLE PRECISION  DIAG(N), FJAC(LDFJAC,N), FVEC(N), QTF(N), R(LR),
     *                  WA1(N), WA2(N), WA3(N), WA4(N), X(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  ACTRED, DELTA, EPSMCH, FNORM, FNORM1, PNORM,
     *                  PRERED, RATIO, SUM, TEMP, XNORM
      INTEGER           I, IFLAG, ITER, J, JM1, L, NCFAIL, NCSUC,
     *                  NSLOW1, NSLOW2
      LOGICAL           JEVAL, SING
C     .. External Functions ..
      DOUBLE PRECISION  F06EJF, X02AJF
      EXTERNAL          F06EJF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          C05NCU, C05NCW, C05NCX, C05NCY, C05NCZ, DGEMV,
     *                  DTPMV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, MOD
C     .. Executable Statements ..
      EPSMCH = X02AJF()
C
      INFO = 1
      IFLAG = 0
      NFEV = 0
      NJEV = 0
C
C     Check the input parameters for errors.
C
      IF (MODE.EQ.2) THEN
         DO 20 J = 1, N
            IF (DIAG(J).LE.ZERO) GO TO 380
   20    CONTINUE
      END IF
C
C     Evaluate the function at the starting point
C     and calculate its norm.
C
      IFLAG = 1
      CALL FCN(N,X,FVEC,FJAC,LDFJAC,IFLAG)
      NFEV = 1
      IF (IFLAG.LT.0) GO TO 380
      FNORM = F06EJF(N,FVEC,1)
C
C     Initialize iteration counter and monitors.
C
      ITER = 1
      NCSUC = 0
      NCFAIL = 0
      NSLOW1 = 0
      NSLOW2 = 0
C
C     Beginning of the outer loop.
C
   40 CONTINUE
      JEVAL = .TRUE.
C
C     Calculate the Jacobian matrix.
C
      IFLAG = 2
      CALL FCN(N,X,FVEC,FJAC,LDFJAC,IFLAG)
      NJEV = NJEV + 1
      IF (IFLAG.LT.0) GO TO 380
C
C     Compute the QR factorization of the Jacobian.
C
      CALL C05NCX(N,FJAC,LDFJAC,WA1,WA2)
C
C     On the first iteration and if MODE is 1, scale according
C     to the norms of the columns of the initial Jacobian.
C
      IF (ITER.EQ.1) THEN
         IF (MODE.NE.2) THEN
            DO 60 J = 1, N
               DIAG(J) = WA2(J)
               IF (WA2(J).EQ.ZERO) DIAG(J) = ONE
   60       CONTINUE
         END IF
C
C        On the first iteration, calculate the norm of the scaled X
C        and initialize the step bound DELTA.
C
         DO 80 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   80    CONTINUE
         XNORM = F06EJF(N,WA3,1)
         DELTA = FACTOR*XNORM
         IF (DELTA.EQ.ZERO) DELTA = FACTOR
      END IF
C
C           T
C     Form Q *FVEC and store in QTF.
C
      DO 100 I = 1, N
         QTF(I) = FVEC(I)
  100 CONTINUE
      DO 160 J = 1, N
         IF (FJAC(J,J).NE.ZERO) THEN
            SUM = ZERO
            DO 120 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
  120       CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 140 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  140       CONTINUE
         END IF
  160 CONTINUE
C
C     Copy the triangular factor of the QR factorization into R.
C
      SING = .FALSE.
      DO 200 J = 1, N
         L = J
         JM1 = J - 1
         IF (J-1.GE.1) THEN
            DO 180 I = 1, JM1
               R(L) = FJAC(I,J)
               L = L + N - I
  180       CONTINUE
         END IF
         R(L) = WA1(J)
         IF (WA1(J).EQ.ZERO) SING = .TRUE.
  200 CONTINUE
C
C     Accumulate the orthogonal factor in FJAC.
C
      CALL C05NCW(N,FJAC,LDFJAC,WA1)
C
C     Rescale if necessary.
C
      IF (MODE.NE.2) THEN
         DO 220 J = 1, N
            DIAG(J) = MAX(DIAG(J),WA2(J))
  220    CONTINUE
      END IF
C
C     Beginning of the inner loop.
C
  240 CONTINUE
C
C     If requested, call FCN to enable printing of iterates.
C
      IF (NPRINT.GT.0) THEN
         IFLAG = 0
         IF (MOD(ITER-1,NPRINT).EQ.0) CALL FCN(N,X,FVEC,FJAC,LDFJAC,
     *       IFLAG)
         IF (IFLAG.LT.0) GO TO 380
      END IF
C
C     Determine the direction P.
C
      CALL C05NCU(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C
C     Store the direction P and X + P. Calculate the norm of P.
C
      DO 260 J = 1, N
         WA1(J) = -WA1(J)
         WA2(J) = X(J) + WA1(J)
         WA3(J) = DIAG(J)*WA1(J)
  260 CONTINUE
      PNORM = F06EJF(N,WA3,1)
C
C     On the first iteration, adjust the initial step bound.
C
      IF (ITER.EQ.1) DELTA = MIN(DELTA,PNORM)
C
C     Evaluate the function at X + P and calculate its norm.
C
      IFLAG = 1
      CALL FCN(N,WA2,WA4,FJAC,LDFJAC,IFLAG)
      NFEV = NFEV + 1
      IF (IFLAG.LT.0) GO TO 380
      FNORM1 = F06EJF(N,WA4,1)
C
C     Compute the scaled actual reduction.
C
      ACTRED = -ONE
      IF (FNORM1.LT.FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C     Compute the scaled predicted reduction.
C
      DO 280 I = 1, N
         WA3(I) = WA1(I)
  280 CONTINUE
      CALL DTPMV('Lower triangle','Transpose','Non-unit diagonal',N,R,
     *           WA3,1)
      DO 300 I = 1, N
         WA3(I) = WA3(I) + QTF(I)
  300 CONTINUE
      TEMP = F06EJF(N,WA3,1)
      PRERED = 1.0D0
      IF (TEMP.LT.FNORM) PRERED = ONE - (TEMP/FNORM)**2
C
C     Compute the ratio of the actual to the predicted
C     reduction.
C
      RATIO = ZERO
      IF (PRERED.GT.ZERO) RATIO = ACTRED/PRERED
C
C     Update the step bound.
C
      IF (RATIO.LT.P1) THEN
         NCSUC = 0
         NCFAIL = NCFAIL + 1
         DELTA = P5*DELTA
      ELSE
         NCFAIL = 0
         NCSUC = NCSUC + 1
         IF (RATIO.GE.P5 .OR. NCSUC.GT.1) DELTA = MAX(DELTA,PNORM/P5)
         IF (ABS(RATIO-ONE).LE.P1) DELTA = PNORM/P5
      END IF
C
C     Test for successful iteration.
C
      IF (RATIO.GE.P0001) THEN
C
C        Successful iteration. Update X, FVEC, and their norms.
C
         DO 320 J = 1, N
            X(J) = WA2(J)
            WA2(J) = DIAG(J)*X(J)
            FVEC(J) = WA4(J)
  320    CONTINUE
         XNORM = F06EJF(N,WA2,1)
         FNORM = FNORM1
         ITER = ITER + 1
      END IF
C
C     Determine the progress of the iteration.
C
      NSLOW1 = NSLOW1 + 1
      IF (ACTRED.GE.P001) NSLOW1 = 0
      IF (JEVAL) NSLOW2 = NSLOW2 + 1
      IF (ACTRED.GE.P1) NSLOW2 = 0
C
C     Test for convergence.
C
      IF (DELTA.LE.XTOL*XNORM .OR. FNORM.EQ.ZERO .OR. PNORM.EQ.ZERO)
     *    INFO = 0
      IF (INFO.NE.1) GO TO 380
C
C     Tests for termination and stringent tolerances.
C
      IF (NFEV.GE.MAXFEV) INFO = 2
      IF (DELTA.LE.EPSMCH*XNORM) INFO = 3
      IF (NSLOW2.EQ.5) INFO = 4
      IF (NSLOW1.EQ.10) INFO = 5
      IF (INFO.NE.1) GO TO 380
C
C     Criterion for recalculating Jacobian.
C
      IF (NCFAIL.NE.2) THEN
C
C        Calculate the rank one modification to the Jacobian
C        and update QTF if necessary.
C
         CALL DGEMV('Transpose',N,N,ONE,FJAC,LDFJAC,WA4,1,ZERO,WA2,1)
         IF (RATIO.GE.P0001) THEN
            DO 340 J = 1, N
               QTF(J) = WA2(J)
  340       CONTINUE
         END IF
         DO 360 J = 1, N
            WA2(J) = (WA2(J)-WA3(J))/PNORM
            WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
  360    CONTINUE
C
C        Compute the QR factorization of the updated Jacobian.
C
         CALL C05NCZ(N,N,R,LR,WA1,WA2,WA3,SING)
         CALL C05NCY(N,N,FJAC,LDFJAC,WA2,WA3)
         CALL C05NCY(1,N,QTF,1,WA2,WA3)
C
C        End of the inner loop.
C
         JEVAL = .FALSE.
         GO TO 240
      END IF
C
C     End of the outer loop.
C
      GO TO 40
  380 CONTINUE
C
C     Termination, either normal or user-imposed.
C
      IF (IFLAG.LT.0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT.GT.0) CALL FCN(N,X,FVEC,FJAC,LDFJAC,IFLAG)
      RETURN
      END
