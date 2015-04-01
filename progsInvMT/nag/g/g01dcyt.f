      DOUBLE PRECISION FUNCTION G01DCY(DXR,D2XR,D3XR,D4XR,D5XR,PR,QR)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     AS APPL. STATIST. ALGORITHM 128.2 (1978), VOL 27
C
C     COMPUTES DAVID-JOHNSON APPROXIMATION FOR THE VARIANCE OF
C     THE RTH LARGEST ORDER STATISTIC FROM THE NORMAL DIST. FOR
C     A SAMPLE SIZE N.
C
C     ARGUMENTS :
C                 DXR - FIRST DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                       EVALUATED AT XR.
C                  :                   :                 :
C                D5XR - FIFTH DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                       EVALUATED AT XR.
C                  PR - EXPECTED VALUE OF RTH LARGEST ORDER STATISTIC
C                       FROM UNIFORM DIST. ( = R/(N+1) R=1,...,N).
C                  QR - 1-PR
C     N.B. XR IS THE INVERSE NORMAL PROBABILITY INTEGRAL OF PR.
C
C     (ANP/AJS)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 D2XR, D3XR, D4XR, D5XR, DXR, PR,
     *                                 QR
C     .. Scalars in Common ..
      DOUBLE PRECISION                 RN2, RN22, RN23
C     .. Local Scalars ..
      DOUBLE PRECISION                 D2XR2, DXR2, FIVETH, FOURTH,
     *                                 HALF, ONEPT5, PRQR, QRMPR, THREE,
     *                                 TWO
C     .. Common blocks ..
      COMMON                           /AG01DC/RN2, RN22, RN23
C     .. Data statements ..
      DATA                             FOURTH/0.25D0/, HALF/0.5D0/,
     *                                 ONEPT5/1.5D0/, FIVETH/
     *                                 1.6666666667D0/, TWO/2.0D0/,
     *                                 THREE/3.0D0/
C     .. Executable Statements ..
      DXR2 = DXR*DXR
      PRQR = PR*QR
      G01DCY = PRQR*DXR2/RN2
C
C     TO ORDER (N+2)**(-2)
C
      QRMPR = QR - PR
      D2XR2 = D2XR*D2XR
      G01DCY = G01DCY + PRQR/RN22*(TWO*QRMPR*DXR*D2XR+PRQR*
     *         (DXR*D3XR+HALF*D2XR2))
C
C     TO ORDER (N+2)**(-3)
C
      G01DCY = G01DCY + PRQR/RN23*(-TWO*QRMPR*DXR*D2XR+
     *         (QRMPR*QRMPR-PRQR)*(TWO*DXR*D3XR+ONEPT5*D2XR2)
     *         +PRQR*QRMPR*(FIVETH*DXR*D4XR+THREE*D2XR*D3XR)
     *         +FOURTH*PRQR*PRQR*(DXR*D5XR+TWO*D2XR*D4XR+FIVETH*D3XR*
     *         D3XR))
      RETURN
      END
