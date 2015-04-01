      SUBROUTINE G01DAX(I,N,CRLN,RINC,SCORE,ERR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes Normal Scores
C
C     A single iteration version of G01DAF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CRLN, ERR, RINC, SCORE
      INTEGER           I, N
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, A3, A4, A5, CK1, CK1A, CK2, CK2A, EARG,
     *                  ENEG, EXP1LN, EXP2LN, F1, F2, PROB1, RIMIN1,
     *                  RLCK1, RLCK2, RLMV1, RLMV2, RLV1, RLV2, RNMINI,
     *                  S, SMALL, T1, T2, T3P, V1, V2, WMV2, X, X1, X2,
     *                  XMIN
      LOGICAL           LOOP1, LOOP2
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02AMF
      EXTERNAL          X01AAF, X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, DBLE, SQRT
C     .. Executable Statements ..
C
C     set up initial values
C
      XMIN = X02AMF()
      F2 = 4.0D0*RINC/3.0D0
      F1 = F2 + F2
      ENEG = LOG(XMIN)
      T1 = ENEG - LOG(F1)
      T2 = ENEG - LOG(F2)
      SMALL = MAX(LOG(X02AJF())*1.5D0,ENEG)
      LOOP1 = .TRUE.
      LOOP2 = .TRUE.
      SCORE = 0.0D0
      T3P = RINC/(3.0D0*SQRT(2.0D0*X01AAF(X)))
      X = 0.0D0
      S = 0.5D0
      A1 = 1.0D0
      EARG = CRLN + DBLE(N-1)*LOG(0.5D0)
      PROB1 = 0.0D0
      IF (EARG.GE.T2) PROB1 = F2*EXP(EARG)
C
C     start of inner loop
C
   20 CONTINUE
      X = X + RINC
      EARG = -0.5D0*X*X
      A2 = 0.0D0
      IF (EARG.GE.ENEG) A2 = EXP(EARG)
      X = X + RINC
      EARG = -0.5D0*X*X
      A3 = 0.0D0
      IF (EARG.GE.ENEG) A3 = EXP(EARG)
      X1 = X
      X = X + RINC
      EARG = -0.5D0*X*X
      A4 = 0.0D0
      IF (EARG.GE.ENEG) A4 = EXP(EARG)
      X = X + RINC
      EARG = -0.5D0*X*X
      A5 = 0.0D0
      IF (EARG.GE.ENEG) A5 = EXP(EARG)
      X2 = X
      V1 = S + (A1+4.0D0*A2+A3)*T3P
      V2 = V1 + (A3+4.0D0*A4+A5)*T3P
      S = V2
      A1 = A5
      WMV2 = 1.0D0 - V2
      IF (V2.LT.1.0D0 .AND. WMV2.GT.0.0D0) THEN
         RLV1 = LOG(V1)
         RLV2 = LOG(V2)
         RLMV1 = LOG(1.0D0-V1)
         RLMV2 = LOG(WMV2)
         EXP1LN = -0.5D0*X1*X1
         EXP2LN = -0.5D0*X2*X2
         RIMIN1 = I - 1
         RNMINI = N - I
         IF (LOOP1) THEN
            RLCK1 = EXP1LN + CRLN + RIMIN1*RLV1 + RNMINI*RLMV1
            RLCK2 = EXP2LN + CRLN + RIMIN1*RLV2 + RNMINI*RLMV2
            IF (RLCK1.GE.SMALL) THEN
               CK1 = 0.0D0
               IF (RLCK1.GE.T1) CK1 = EXP(RLCK1)*F1
               CK2 = 0.0D0
               IF (RLCK2.GE.T2) CK2 = EXP(RLCK2)*F2
               SCORE = SCORE + CK1*X1 + CK2*X2
               PROB1 = PROB1 + CK1 + CK2
            ELSE
               LOOP1 = .FALSE.
            END IF
         END IF
         IF (LOOP2) THEN
            RLCK1 = EXP1LN + CRLN + RIMIN1*RLMV1 + RNMINI*RLV1
            RLCK2 = EXP2LN + CRLN + RIMIN1*RLMV2 + RNMINI*RLV2
            IF (RLCK2.GE.RLCK1 .OR. RLCK1.GE.SMALL) THEN
               CK1 = 0.0D0
               IF (RLCK1.GE.T1) CK1 = EXP(RLCK1)*F1
               CK2 = 0.0D0
               IF (RLCK2.GE.T2) CK2 = EXP(RLCK2)*F2
               CK1A = 0.0D0
               IF (X1.GE.1.0D0) THEN
                  CK1A = CK1*X1
               ELSE IF (CK1.GE.XMIN/X1) THEN
                  CK1A = CK1*X1
               END IF
               CK2A = 0.0D0
               IF (X2.GE.1.0D0) THEN
                  CK2A = CK2*X2
               ELSE IF (CK2.GE.XMIN/X2) THEN
                  CK2A = CK2*X2
               END IF
               SCORE = SCORE - CK1A - CK2A
               PROB1 = PROB1 + CK1 + CK2
            ELSE
               LOOP2 = .FALSE.
            END IF
         END IF
         IF (LOOP1 .OR. LOOP2) GO TO 20
      END IF
      ERR = ABS(1.0D0-PROB1)*X
      RETURN
      END
