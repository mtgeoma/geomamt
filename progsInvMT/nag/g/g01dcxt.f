      DOUBLE PRECISION FUNCTION G01DCX(DXR,D2XR,D3XR,D4XR,D5XR,PR,QR,
     *                                 DXS,D2XS,D3XS,D4XS,D5XS,PS)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     AS APPL. STATIST. ALGORITHM AS 128.3 (1978), VOL. 27
C     DAVIS C.S. AND STEPHENS M.A.
C
C     COMPUTES DAVID-JOHNSON APPROXIMATION FOR COVARIANCE BETWEEN RTH
C     AND STH ORDER STATISTICS FROM THE NORMAL DIST. FOR A SAMPLE
C     SIZE N.
C
C     ARGUMENTS :
C                 DXR - FIRST DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                       EVALUATED AT XR.
C                  :                     :               :
C                D5XR - FIFTH DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                       EVALUATED AT XR.
C                  PR - EXPECTED VALUE OF RTH ORDER STATISTIC FROM
C                       UNIFORM DIST. ( = R/(N+1) R=1,...N ).
C                       QR - 1-PR.
C                 DXS - FIRST DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                       EVALUATED AT XS.
C                  :                     :             :
C                D5XS - FIFTH DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                       EVALUATED AT XS.
C                  PS - EXPECTED VALUE OF STH ORDER STATISTIC FROM
C                       UNIFORM DISTRIBUTION. ( = S/(N+1) S=1,...,N).
C
C     N.B. XR IS THE INVERSE NORMAL PROBABILITY INTEGRAL OF PR ETC.
C
C     (ANP/AJS)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 D2XR, D2XS, D3XR, D3XS, D4XR,
     *                                 D4XS, D5XR, D5XS, DXR, DXS, PR,
     *                                 PS, QR
C     .. Scalars in Common ..
      DOUBLE PRECISION                 RN2, RN22, RN23
C     .. Local Scalars ..
      DOUBLE PRECISION                 EIGTH, FIVE6, FOURTH, HALF, ONE,
     *                                 ONEPT5, PR2, PRQR, PRQS, PS2,
     *                                 PSQR, PSQS, QR2, QRMPR, QS, QS2,
     *                                 QSMPS, TERM1, TERM2, TERM3,
     *                                 TERM4, TERM5, THREE, TWELTH, TWO
C     .. Common blocks ..
      COMMON                           /AG01DC/RN2, RN22, RN23
C     .. Data statements ..
      DATA                             TWELTH/0.0833333333333D0/,
     *                                 EIGTH/0.125D0/, FOURTH/0.25D0/,
     *                                 HALF/0.5D0/, FIVE6/
     *                                 0.833333333333D0/, ONE/1.0D0/,
     *                                 ONEPT5/1.5D0/, TWO/2.0D0/,
     *                                 THREE/3.0D0/
C     .. Executable Statements ..
      QS = ONE - PS
      PRQS = PR*QS
      G01DCX = PRQS*DXR*DXS/RN2
C
C     TO ORDER (N+2)**(-2)
C
      QRMPR = QR - PR
      QSMPS = QS - PS
      PRQR = PR*QR
      PSQS = PS*QS
      G01DCX = G01DCX + PRQS/RN22*(QRMPR*D2XR*DXS+QSMPS*DXR*D2XS+HALF*
     *         PRQR*D3XR*DXS+HALF*PSQS*DXR*D3XS+HALF*PRQS*D2XR*D2XS)
C
C     TO ORDER (N+2)**(-3)
C
      PR2 = PR*PR
      QR2 = QR*QR
      PS2 = PS*PS
      QS2 = QS*QS
      PSQR = PS*QR
      TERM1 = -D2XR*DXS*QRMPR - QSMPS*DXR*D2XS + (QRMPR*QRMPR-PRQR)
     *        *D3XR*DXS
      TERM2 = (QSMPS*QSMPS-PSQS)*DXR*D3XS +
     *        (ONEPT5*QRMPR*QSMPS+HALF*PSQR-TWO*PRQS)*D2XR*D2XS
      TERM3 = FIVE6*(PRQR*QRMPR*D4XR*DXS+PSQS*QSMPS*DXR*D4XS) +
     *        (PRQS*QRMPR+HALF*PRQR*QSMPS)*D3XR*D2XS
      TERM4 = (PRQS*QSMPS+HALF*PSQS*QRMPR)*D2XR*D3XS +
     *        EIGTH*(PR2*QR2*D5XR*DXS+PS2*QS2*DXR*D5XS)
      TERM5 = FOURTH*(PR2*QR*QS*D4XR*D2XS+PR*PS*QS2*D2XR*D4XS) +
     *        TWELTH*(TWO*PR2*QS2+THREE*PR*QR*PS*QS)*D3XR*D3XS
      G01DCX = G01DCX + PRQS/RN23*(TERM1+TERM2+TERM3+TERM4+TERM5)
      RETURN
      END
