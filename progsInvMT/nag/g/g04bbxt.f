      SUBROUTINE G04BBX(N,IBLOCK,KBLOCK,NT,IT,C,LDC,IREP,EF,NTDF,TMEAN,
     +                  ACC,WK,IWARN)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     MARK 17 REVISED. IER-1663 (JUN 1995).
C
C     Computes treatment effects for non-orthogonal designs and the
C     generalised inverse for the reduced treatment matrix
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION ACC
      INTEGER IBLOCK,IWARN,KBLOCK,LDC,N,NT,NTDF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION C(LDC,NT),EF(NT),TMEAN(NT),WK(3*NT)
      INTEGER IREP(NT),IT(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACCE,RMEAN,VAR
      INTEGER I,ICOL,IFAULT,II,INCI,INCJ,J,JJ,JTH,K,KK,KTH,NBLOCK
C     ..
C     .. External Functions ..
      DOUBLE PRECISION F06EAF
      EXTERNAL F06EAF
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,F02FAF,F06PAF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,SQRT
C     ..
      IWARN = 0
      DO 40 I = 1,NT
          DO 20 J = 1,NT
              C(J,I) = 0.0D0
   20     CONTINUE
   40 CONTINUE
C
C      Compute the NN' matrix
C
      NBLOCK = ABS(IBLOCK)
      IF (IBLOCK.LT.0) THEN
          INCI = 1
          INCJ = NBLOCK
      ELSE
          INCI = KBLOCK
          INCJ = 1
      END IF
      II = 1
      DO 120 I = 1,NBLOCK
          JJ = II
          DO 100 J = 1,KBLOCK
              JTH = IT(JJ)
              KK = II
              DO 60 K = 1,J - 1
                  IF (JTH.EQ.IT(KK)) THEN
                      C(JTH,JTH) = C(JTH,JTH) - 1.0D0
                  END IF
                  KK = KK + INCJ
   60         CONTINUE
              C(JTH,JTH) = C(JTH,JTH) - 1.0D0
              KK = KK + INCJ
              DO 80 K = J + 1,KBLOCK
                  KTH = IT(KK)
                  IF (JTH.LE.KTH) THEN
                      C(KTH,JTH) = C(KTH,JTH) - 1.0D0
                  ELSE
                      C(JTH,KTH) = C(JTH,KTH) - 1.0D0
                  END IF
                  KK = KK + INCJ
   80         CONTINUE
              JJ = JJ + INCJ
  100     CONTINUE
          II = II + INCI
  120 CONTINUE
C
C     Compute A matrix
C
      DO 160 I = 1,NT
          C(I,I) = IREP(I) + C(I,I)/DBLE(KBLOCK)
          DO 140 J = I + 1,NT
              C(J,I) = C(J,I)/DBLE(KBLOCK)
  140     CONTINUE
  160 CONTINUE
C
C     Find efficiency factors and parameter estimates
C
C     Compute eigenvalues of A
C
      IFAULT = 1
      CALL F02FAF('V','L',NT,C,LDC,EF,WK,3*NT,IFAULT)
      IF (IFAULT.NE.0) THEN
          IWARN = -4
          GO TO 320
      END IF
C
C     Check for zero values and compute q*U*INV(E)
C
      NTDF = 0
      ACCE = ACC*EF(NT)
      DO 180 I = 1,NT
          IF (EF(I).GT.ACCE) THEN
              NTDF = NTDF + 1
              WK(I) = F06EAF(NT,C(1,I),1,TMEAN,1)/EF(I)
          ELSE
              EF(I) = 0.0D0
          END IF
  180 CONTINUE
      IF (NTDF.EQ.0) THEN
          IWARN = -3
          GO TO 320
      END IF
      ICOL = NT - NTDF + 1
      CALL F06PAF('N',NT,NTDF,1.0D0,C(1,ICOL),LDC,WK(ICOL),1,0.0D0,
     +            TMEAN,1)
C
C     Compute generalised inverse
C
      DO 220 I = NT,1,-1
          DO 200 J = ICOL,NT
              WK(J) = C(I,J)/EF(J)
  200     CONTINUE
          CALL F06PAF('N',I,NTDF,1.0D0,C(1,ICOL),LDC,WK(ICOL),1,0.0D0,
     +                WK(NT+1),1)
          CALL DCOPY(I,WK(NT+1),1,C(I,1),LDC)
  220 CONTINUE
      DO 230 I = 1,NT - 1
          CALL DCOPY(NT-I,C(I+1,I),1,C(I,I+1),LDC)
  230 CONTINUE
C
C     Compute se of differences in means
C
      DO 260 I = 1,NT - 1
          DO 240 J = I + 1,NT
              VAR = C(I,I) + C(J,J) - 2.0D0*C(I,J)
              IF (VAR.GT.0.0D0) THEN
                  C(J,I) = SQRT(VAR)
              ELSE
                  C(J,I) = 0.0D0
                  IWARN = -1
              END IF
  240     CONTINUE
  260 CONTINUE
C
C     Scale efficiency factors
C
      RMEAN = 0.0D0
      DO 280 I = 1,NT
          RMEAN = RMEAN + IREP(I)
  280 CONTINUE
      RMEAN = RMEAN/DBLE(NT)
      DO 300 I = ICOL,NT
          EF(I) = EF(I)/RMEAN
  300 CONTINUE
      IF (NTDF.LT.NT-1) IWARN = -2
  320 CONTINUE
      RETURN
      END
