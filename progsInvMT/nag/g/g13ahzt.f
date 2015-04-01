      SUBROUTINE G13AHZ(ST,NST,NP,ND,NQ,NPS,NDS,NQS,NS,PHI,THETA,SPHI,
     *                  STHETA,NPAR,C,RMS,NFV,FVA,FSD,AEX,AAL,AEXR)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AHZ CARRIES OUT THE CALCULATIONS INVOLVED
C     IN FORECASTING
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, RMS
      INTEGER           ND, NDS, NFV, NP, NPAR, NPS, NQ, NQS, NS, NST
C     .. Array Arguments ..
      DOUBLE PRECISION  AAL(NST), AEX(NST), AEXR(NST), FSD(NFV),
     *                  FVA(NFV), PHI(NPAR), SPHI(NPAR), ST(NST),
     *                  STHETA(NPAR), THETA(NPAR)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, AL, D, U, UX, X, ZERO
      INTEGER           I, J, K, KFA, LDS, NSTM, NT
C     .. External Subroutines ..
      EXTERNAL          G13AHX, G13AHY
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, U/1.0D0/
C     .. Executable Statements ..
C
C     CALL THE SUBROUTINE WHICH DECOMPOSES THE STATE SET INTO
C     THREE CONSTITUENT ARRAYS  AEX, AAL AND AEXR
C
      CALL G13AHY(ST,NST,NP,ND,NQ,NPS,NDS,NQS,NS,AEX,AAL,AEXR)
C
C     DEFINE POINT OF SUBDIVISION OF AEX, AAL AND AEXR
C
      LDS = NST - ND - NDS*NS
      NT = LDS + 1
      NSTM = NST - 1
      D = C
      A = ZERO
      KFA = 0
C
C     PROCESS ONE FORECAST VALUE AT A TIME
C
   20 DO 280 I = 1, NFV
C
C        GO THROUGH THE CALCULATIONS WHICH OBTAIN THE FORECAST
C        OF THE NEXT DIFFERENCED VALUE
C
         AL = A
         IF (NP.LE.0) GO TO 60
         DO 40 J = 1, NP
            K = NT - J
            AL = AL + PHI(J)*AAL(K)
   40    CONTINUE
   60    IF (NQ.LE.0) GO TO 100
         DO 80 J = 1, NQ
            K = NT - J
            AL = AL - THETA(J)*AEXR(K)
   80    CONTINUE
  100    X = AL + C
         IF (NPS.LE.0) GO TO 140
         DO 120 J = 1, NPS
            K = NT - J*NS
            X = X + SPHI(J)*(AEX(K)-C)
  120    CONTINUE
  140    IF (NQS.LE.0) GO TO 180
         DO 160 J = 1, NQS
            K = NT - J*NS
            X = X - STHETA(J)*AAL(K)
  160    CONTINUE
C
C        CALL THE SUBROUTINE WHICH OBTAINS THE UNDIFFERENCED VALUE
C        RELATING TO THE FORECAST DIFFERENCED VALUE
C
  180    CALL G13AHX(AEX,NST,ND,NDS,NS,1,UX,X)
C
C        STORE UNDIFFERENCED VALUE AS A FORECAST VALUE OR AS ITS S.D.
C
         IF (KFA.EQ.0) GO TO 200
         FSD(I) = UX
         GO TO 220
  200    FVA(I) = UX
  220    IF (I.EQ.NFV) GO TO 280
C
C        SHIFT THE AEX, AAL AND AEXR ARRAYS UP ONE VALUE, AND PUT
C        NEW VALUES AT END OF AEX, AND AT SUBDIVISION POINT IN AAL AND
C        AEXR
C
         IF (NSTM.LE.0) GO TO 260
         DO 240 J = 1, NSTM
            AEX(J) = AEX(J+1)
            AAL(J) = AAL(J+1)
            AEXR(J) = AEXR(J+1)
  240    CONTINUE
  260    AEX(NST) = UX
         AAL(LDS) = AL
         AEXR(LDS) = A
         A = ZERO
  280 CONTINUE
C
C     REPEAT LOOP WITH PARAMETERS WHICH WILL GIVE S.D. INSTEAD
C     OF FORECAST VALUES
C
      IF (KFA.EQ.1) GO TO 320
      A = U
      KFA = 1
      C = ZERO
      DO 300 J = 1, NST
         AEX(J) = ZERO
         AAL(J) = ZERO
         AEXR(J) = ZERO
  300 CONTINUE
      GO TO 20
  320 C = D
      D = ZERO
      DO 340 I = 1, NFV
         D = D + FSD(I)**2
         FSD(I) = SQRT(RMS*D)
  340 CONTINUE
      RETURN
      END
