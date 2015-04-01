      DOUBLE PRECISION FUNCTION G01DCV(N)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATES AN APPROXIMATION TO THE VARIANCE OF
C     THE LARGEST NORMAL ORDER STATISTIC IN A SAMPLE OF
C     SIZE N.
C
C     ARGUMENTS :
C                 N - SAMPLE SIZE
C
C     (ANP/AJS)
C
C     .. Scalar Arguments ..
      INTEGER                          N
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, AN, AN1, AN2, AN3, AN4, B1,
     *                                 B2, B3, B4, B5, ETA
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG
C     .. Data statements ..
      DATA                             B1/0.2737885D-1/,
     *                                 B2/-0.62877655D0/,
     *                                 B3/0.81480712D-1/,
     *                                 B4/-0.6835144D-2/,
     *                                 B5/0.2607097D-3/, A/0.45D-1/
C     .. Executable Statements ..
      AN = A + N
      AN1 = LOG(AN)
      AN2 = AN1*AN1
      AN3 = AN1*AN2
      AN4 = AN1*AN3
      ETA = B1 + B2*AN1 + B3*AN2 + B4*AN3 + B5*AN4
      G01DCV = EXP(ETA)
      RETURN
      END
