      DOUBLE PRECISION FUNCTION G01DBZ(I,N)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     AS APPL. STATIST. ALGORITHM 177.4 (1982), VOL. 31
C
C     CALCULATES CORRECTION FOR TAIL AREA OF NORMAL DISTRIBUTION
C     CORRESPONDING TO ITH LARGEST RANKIT IN SAMPLE SIZE N.
C
C     ARGUMENTS : I - ITH LARGEST ORDER STATISTIC (FOR I.LE.7)
C                 N - SAMPLE SIZE (FOR N.LE.20)
C
C     .. Scalar Arguments ..
      INTEGER                          I, N
C     .. Local Scalars ..
      DOUBLE PRECISION                 AMIC, AN, C14, ONE, ZERO
C     .. Local Arrays ..
      DOUBLE PRECISION                 C1(7), C2(7), C3(7)
C     .. Data statements ..
      DATA                             C1(1), C1(2), C1(3), C1(4),
     *                                 C1(5), C1(6), C1(7)/9.5D0,
     *                                 28.7D0, 1.9D0, 0.0D0, -7.0D0,
     *                                 -6.2D0, -1.6D0/, C2(1), C2(2),
     *                                 C2(3), C2(4), C2(5), C2(6),
     *                                 C2(7)/-6.195D3, -9.569D3,
     *                                 -6.728D3, -17.614D3, -8.278D3,
     *                                 -3.570D3, 1.075D3/, C3(1), C3(2),
     *                                 C3(3), C3(4), C3(5), C3(6),
     *                                 C3(7)/9.338D4, 1.7516D5,
     *                                 4.1040D5, 2.157D6, 2.376D6,
     *                                 2.065D6, 2.065D6/, AMIC/1.0D-6/,
     *                                 C14/1.9D-5/, ZERO/0.0D0/,
     *                                 ONE/1.0D0/
C     .. Executable Statements ..
      G01DBZ = C14
      IF (I*N.EQ.4) RETURN
      G01DBZ = ZERO
      IF (I.GT.7) RETURN
      IF (I.NE.4 .AND. N.GT.20) RETURN
      IF (I.EQ.4 .AND. N.GT.40) RETURN
      AN = N
      AN = ONE/(AN*AN)
      G01DBZ = (C1(I)+AN*(C2(I)+AN*C3(I)))*AMIC
      RETURN
      END
