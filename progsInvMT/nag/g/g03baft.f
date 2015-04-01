      SUBROUTINE G03BAF(STAND,G,NVAR,K,FL,LDF,FLR,R,LDR,ACC,MAXIT,ITER,
     *                  WK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1036 (JUN 1993).
C
C     Performs orthogonal rotations.
C
C     The NVAR x K matrix in FL is rotated to give matrix in FLR using
C     the rotations in R.
C     If G = 1 varimax rotation is used
C     If G = 0 quartimax rotation is used
C     If STAND = 'S' or 's' the rotations are computed from the
C     row-standardized matrix.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC, G
      INTEGER           IFAIL, ITER, K, LDF, LDR, MAXIT, NVAR
      CHARACTER         STAND
C     .. Array Arguments ..
      DOUBLE PRECISION  FL(LDF,K), FLR(LDF,K), R(LDR,K),
     *                  WK(2*NVAR+K*K+5*(K-1))
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, B1, C, C1, CP, D1, EPS, F2, FI, FJ, OSQRT2,
     *                  Q, QD, QN, QNS, QP, S, SCALE, SP, SUM, TEMP, U,
     *                  V
      INTEGER           I, IERROR, IFAULT, J, L, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP1(1), WKSP2(1)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DNRM2, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, F02WEF, F06BAF,
     *                  F06QFF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, SIGN, SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (K.LT.2) THEN
         WRITE (P01REC,FMT=99993) K
      ELSE IF (NVAR.LT.K) THEN
         WRITE (P01REC,FMT=99992) NVAR, K
      ELSE IF (LDF.LT.NVAR) THEN
         WRITE (P01REC,FMT=99999) LDF, NVAR
      ELSE IF (LDR.LT.K) THEN
         WRITE (P01REC,FMT=99998) LDR, K
      ELSE IF (STAND.NE.'U' .AND. STAND.NE.'S' .AND. STAND.NE.'u' .AND.
     *         STAND.NE.'s') THEN
         WRITE (P01REC,FMT=99995) STAND
      ELSE IF (G.LT.0.0D0) THEN
         WRITE (P01REC,FMT=99994) G
      ELSE IF (ACC.LT.0.0D0) THEN
         WRITE (P01REC,FMT=99997) ACC
      ELSE IF (MAXIT.LE.0) THEN
         WRITE (P01REC,FMT=99990) MAXIT
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         OSQRT2 = 1.0D0/SQRT(2.0D0)
         EPS = X02AJF()
         IF (ACC.GT.EPS .AND. ACC.LT.1.0D0) EPS = ACC
C
C        STANDARDISE BY ROWS IF STAND=S
C
         IF (STAND.EQ.'S' .OR. STAND.EQ.'s') THEN
            DO 20 I = 1, NVAR
               WK(I) = DNRM2(K,FL(I,1),LDF)
               SCALE = 1.0D0/WK(I)
               CALL DSCAL(K,SCALE,FL(I,1),LDF)
   20       CONTINUE
         END IF
         ITER = 0
         QP = 0.0D0
         SCALE = -G/DBLE(NVAR)
         CALL F06QFF('G',NVAR,K,FL,LDF,FLR,LDF)
C
C        Start iterations
C
   40    ITER = ITER + 1
C
C        Use Cooley-lohnes algorithm for initial estimates
C
         Q = 0.0D0
         DO 80 J = 1, K
            FJ = 0.0D0
            DO 60 I = 1, NVAR
               FI = FLR(I,J)*FLR(I,J)
               Q = Q + FI*FI
               FJ = FJ + FI
   60       CONTINUE
            Q = Q - SCALE*FJ
   80    CONTINUE
         IF (ABS(Q-QP).GT.EPS*MAX(QP,1.0D0)) THEN
            QP = ABS(Q)
            DO 160 I = 1, K - 1
               DO 140 J = I + 1, K
                  A1 = 0.0D0
                  B1 = 0.0D0
                  C1 = 0.0D0
                  D1 = 0.0D0
                  DO 100 L = 1, NVAR
                     FI = FLR(L,I)
                     FJ = FLR(L,J)
                     U = FI*FI - FJ*FJ
                     V = 2.0D0*FI*FJ
                     A1 = A1 + U
                     B1 = B1 + V
                     C1 = C1 + U*U - V*V
                     D1 = D1 + 2.0D0*U*V
  100             CONTINUE
                  QN = D1 + 2.0D0*A1*B1*SCALE
                  QD = C1 + (A1*A1-B1*B1)*SCALE
                  QNS = SIGN(1.0D0,QN)
                  IF (QD.GE.0.0D0) THEN
                     CALL F06BAF(QD,QN,C,S)
                     TEMP = SQRT((1.0D0+C)/2.0D0)
                     CP = SQRT((1.0D0+TEMP)/2.0D0)
                     SP = S/(4.0D0*CP*TEMP)
                  ELSE
                     CALL F06BAF(QD,QN,C,S)
                     TEMP = SQRT((1.0D0+C)/2.0D0)
                     C = SQRT((1.0D0+TEMP)/2.0D0)
                     S = S/(4.0D0*C*TEMP)
                     IF (QNS.GT.0) THEN
                        SP = OSQRT2*(C+S)
                        CP = OSQRT2*(C-S)
                     ELSE
                        CP = OSQRT2*(C+S)
                        SP = -OSQRT2*(C-S)
                     END IF
                  END IF
                  DO 120 L = 1, NVAR
                     TEMP = FLR(L,I)*CP + FLR(L,J)*SP
                     FLR(L,J) = FLR(L,J)*CP - FLR(L,I)*SP
                     FLR(L,I) = TEMP
  120             CONTINUE
  140          CONTINUE
  160       CONTINUE
            IF (ITER.LE.MAXIT) GO TO 40
            IERROR = 3
            WRITE (P01REC,FMT=99996)
         END IF
C
C        Use Lawley-Maxwell method for final improvement and
C        calculation of rotation matrix.
C
         DO 200 J = 1, K
            SUM = 0.0D0
            DO 180 I = 1, NVAR
               F2 = FLR(I,J)*FLR(I,J)
               WK(NVAR+I) = FLR(I,J)*F2
               SUM = SUM + F2
  180       CONTINUE
            TEMP = SUM*SCALE
            CALL DAXPY(NVAR,TEMP,FLR(1,J),1,WK(NVAR+1),1)
            CALL DGEMV('T',NVAR,K,1.0D0,FL,LDF,WK(NVAR+1),1,0.0D0,
     *                 FLR(1,J),1)
  200    CONTINUE
         IFAULT = 1
         CALL F02WEF(K,K,FLR,LDF,0,WKSP1,1,.TRUE.,WKSP2,1,WK(NVAR+1),
     *               .TRUE.,R,LDR,WK(2*NVAR+1),IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 2
            WRITE (P01REC,FMT=99991)
         ELSE
            DO 220 J = 1, K
               CALL DCOPY(K,R(1,J),1,WK(NVAR+1),1)
               CALL DGEMV('N',K,K,1.0D0,FLR,LDF,WK(NVAR+1),1,0.0D0,
     *                    R(1,J),1)
  220       CONTINUE
            DO 240 J = 1, K
               CALL DGEMV('N',NVAR,K,1.0D0,FL,LDF,R(1,J),1,0.0D0,
     *                    FLR(1,J),1)
  240       CONTINUE
         END IF
C
C        De-Normalize rows of FLR
C
         IF (STAND.EQ.'S' .OR. STAND.EQ.'s') THEN
            DO 260 I = 1, NVAR
               CALL DSCAL(K,WK(I),FL(I,1),LDF)
               CALL DSCAL(K,WK(I),FLR(I,1),LDF)
  260       CONTINUE
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, LDF.lt.NVAR: LDF = ',I16,' NVAR = ',I16)
99998 FORMAT (' ** On entry, LDR.lt.K: LDR = ',I16,' K = ',I16)
99997 FORMAT (' ** On entry, ACC.lt.0.0: ACC = ',D13.5)
99996 FORMAT (' ** Procedure failed to converge in MAXIT iterations: M',
     *       'AXIT =',I16)
99995 FORMAT (' ** On entry, STAND is not a valid character: STAND = ',
     *       A1)
99994 FORMAT (' ** On entry, G.lt.0.0: G = ',D13.5)
99993 FORMAT (' ** On entry, K.lt.2: K = ',I16)
99992 FORMAT (' ** On entry, NVAR.lt.K: NVAR = ',I16,' K = ',I16)
99991 FORMAT (' ** SVD failed to converge')
99990 FORMAT (' ** On entry, MAXIT.le.0: MAXIT = ',I16)
      END
