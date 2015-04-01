      SUBROUTINE C05PDZ(IREVCM,N,X,FVEC,FJAC,LDFJAC,XTOL,DIAG,MODE,
     *                  FACTOR,INFO,R,LR,QTF,WA1,WA2,WA3,WA4)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     **********
C
C     SUBROUTINE C05PDZ (based on MINPACK routine HYBRJ )
C
C     The purpose of C05PDZ is to find a zero of a system of
C     N nonlinear functions in N variables by a modification
C     of the POWELL hybrid method.
C
C     **********
C
C     Revised to call BLAS.
C     P.J.D. Mayes, J.J. Du Croz, NAG Central Office, September 1987.
C     Revised to reverse communication.
C     M.S. Derakhshan, NAG Central Office, September 1988.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, P1, P5, P001, P0001, ZERO
      PARAMETER         (ONE=1.0D0,P1=0.1D0,P5=0.5D0,P001=0.001D0,
     *                  P0001=0.0001D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FACTOR, XTOL
      INTEGER           INFO, IREVCM, LDFJAC, LR, MODE, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DIAG(N), FJAC(LDFJAC,N), FVEC(N), QTF(N), R(LR),
     *                  WA1(N), WA2(N), WA3(N), WA4(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACTRED, DELTA, EPSMCH, FNORM, FNORM1, PNORM,
     *                  PRERED, RATIO, SUM, TEMP, XNORM
      INTEGER           I, IREV1, IREV2, ITER, J, JM1, L, NCFAIL, NCSUC,
     *                  NSLOW1, NSLOW2
      LOGICAL           JEVAL, SING
C     .. External Functions ..
      DOUBLE PRECISION  F06EJF, X02AJF
      EXTERNAL          F06EJF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          C05NCU, C05NCW, C05NCX, C05NCY, C05NCZ, DGEMV,
     *                  DTPMV, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Save statement ..
      SAVE
C     .. Executable Statements ..
      GO TO (300,40,80) IREVCM
      IREV1 = 0
      IREV2 = 0
      EPSMCH = X02AJF()
C
      INFO = 1
C
C     Check the input parameters for errors.
C
      IF (MODE.EQ.2) THEN
         DO 20 J = 1, N
            IF (DIAG(J).LE.ZERO) GO TO 460
   20    CONTINUE
      END IF
C
C     Return to user-level to evaluate the function
C     at the starting point and calculate its norm.
C
      IREVCM = 2
      RETURN
   40 CONTINUE
      IF (IREV2.EQ.1) GO TO 340
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
   60 CONTINUE
      JEVAL = .TRUE.
C
C     Return to user-level to calculate the Jacobian matrix.
C
      IREVCM = 3
      RETURN
   80 CONTINUE
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
            DO 100 J = 1, N
               DIAG(J) = WA2(J)
               IF (WA2(J).EQ.ZERO) DIAG(J) = ONE
  100       CONTINUE
         END IF
C
C        On the first iteration, calculate the norm of the scaled X
C        and initialize the step bound DELTA.
C
         DO 120 J = 1, N
            WA3(J) = DIAG(J)*X(J)
  120    CONTINUE
         XNORM = F06EJF(N,WA3,1)
         DELTA = FACTOR*XNORM
         IF (DELTA.EQ.ZERO) DELTA = FACTOR
      END IF
C
C           T
C     Form Q *FVEC and store in QTF.
C
      DO 140 I = 1, N
         QTF(I) = FVEC(I)
  140 CONTINUE
      DO 200 J = 1, N
         IF (FJAC(J,J).NE.ZERO) THEN
            SUM = ZERO
            DO 160 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
  160       CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 180 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  180       CONTINUE
         END IF
  200 CONTINUE
C
C     Copy the triangular factor of the QR factorization into R.
C
      SING = .FALSE.
      DO 240 J = 1, N
         L = J
         JM1 = J - 1
         IF (J-1.GE.1) THEN
            DO 220 I = 1, JM1
               R(L) = FJAC(I,J)
               L = L + N - I
  220       CONTINUE
         END IF
         R(L) = WA1(J)
         IF (WA1(J).EQ.ZERO) SING = .TRUE.
  240 CONTINUE
C
C     Accumulate the orthogonal factor in FJAC.
C
      CALL C05NCW(N,FJAC,LDFJAC,WA1)
C
C     Rescale if necessary.
C
      IF (MODE.NE.2) THEN
         DO 260 J = 1, N
            DIAG(J) = MAX(DIAG(J),WA2(J))
  260    CONTINUE
      END IF
C
C     Beginning of the inner loop.
C
  280 CONTINUE
C
C     If requested, return to user-level to print the iterates.
C
      IREVCM = 1
      RETURN
  300 CONTINUE
      IF (IREV1.EQ.1) GO TO 480
C
C     Determine the direction P.
C
      CALL C05NCU(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C
C     Store the direction P and X + P. Calculate the norm of P.
C
      DO 320 J = 1, N
         WA1(J) = -WA1(J)
         WA2(J) = X(J) + WA1(J)
         WA3(J) = DIAG(J)*WA1(J)
  320 CONTINUE
      PNORM = F06EJF(N,WA3,1)
C
C     On the first iteration, adjust the initial step bound.
C
      IF (ITER.EQ.1) DELTA = MIN(DELTA,PNORM)
C
C     Return to user-level to evaluate the function
C     at X + P and calculate its norm.
C
      IREVCM = 2
      IREV2 = 1
      CALL DSWAP(N,WA2,1,X,1)
      CALL DSWAP(N,WA4,1,FVEC,1)
      RETURN
  340 CONTINUE
      CALL DSWAP(N,X,1,WA2,1)
      CALL DSWAP(N,FVEC,1,WA4,1)
      FNORM1 = F06EJF(N,WA4,1)
C
C     Compute the scaled actual reduction.
C
      ACTRED = -ONE
      IF (FNORM1.LT.FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C     Compute the scaled predicted reduction.
C
      DO 360 I = 1, N
         WA3(I) = WA1(I)
  360 CONTINUE
      CALL DTPMV('Lower triangle','Transpose','Non-unit diagonal',N,R,
     *            WA3,1)
      DO 380 I = 1, N
         WA3(I) = WA3(I) + QTF(I)
  380 CONTINUE
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
         DO 400 J = 1, N
            X(J) = WA2(J)
            WA2(J) = DIAG(J)*X(J)
            FVEC(J) = WA4(J)
  400    CONTINUE
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
      IF (INFO.NE.1) GO TO 460
C
C     Tests for termination and stringent tolerances.
C
      IF (DELTA.LE.EPSMCH*XNORM) INFO = 3
      IF (NSLOW2.EQ.5) INFO = 4
      IF (NSLOW1.EQ.10) INFO = 5
      IF (INFO.NE.1) GO TO 460
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
            DO 420 J = 1, N
               QTF(J) = WA2(J)
  420       CONTINUE
         END IF
         DO 440 J = 1, N
            WA2(J) = (WA2(J)-WA3(J))/PNORM
            WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
  440    CONTINUE
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
         GO TO 280
      END IF
C
C     End of the outer loop.
C
      GO TO 60
  460 CONTINUE
C
C     Termination, either normal or user-imposed.
C
      IREVCM = 1
      IREV1 = 1
      RETURN
  480 CONTINUE
      IREVCM = 0
      RETURN
      END
