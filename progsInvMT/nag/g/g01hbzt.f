      DOUBLE PRECISION FUNCTION G01HBZ(X)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Calculates the multivariate normal density
C     Called by D01AJF, used when N = 3
C
C     .. Parameters ..
      INTEGER                          NMAX
      PARAMETER                        (NMAX=10)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Scalars in Common ..
      DOUBLE PRECISION                 CONST, L1, L2, RHO, SD1, SD2, U1,
     *                                 U2
      INTEGER                          IND
C     .. Arrays in Common ..
      DOUBLE PRECISION                 C(NMAX,NMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION                 AM1, AM2, E, P, UFLOW, ZL1, ZL2,
     *                                 ZU1, ZU2
      INTEGER                          IFAULT
C     .. External Functions ..
      DOUBLE PRECISION                 G01HAF, X02AMF
      EXTERNAL                         G01HAF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG, MAX, MIN
C     .. Common blocks ..
      COMMON                           /AG01HB/C, CONST, L1, L2, U1, U2,
     *                                 SD1, SD2, RHO, IND
C     .. Executable Statements ..
      UFLOW = LOG(X02AMF())
      E = -0.5D0*X*X - CONST
      IF (E.GT.UFLOW) THEN
C
C        compute probability for conditional bivariate
C
         AM1 = C(2,1)*X
         AM2 = C(3,1)*X
         IF (IND.EQ.0) THEN
            ZL1 = (L1-AM1)/SD1
            ZL2 = (L2-AM2)/SD2
            ZU1 = (U1-AM1)/SD1
            ZU2 = (U2-AM2)/SD2
            IFAULT = 0
            P = G01HAF(ZU1,ZU2,RHO,IFAULT) + G01HAF(ZL1,ZL2,RHO,IFAULT)
     *           - G01HAF(ZL1,ZU2,RHO,IFAULT) - G01HAF(ZU1,ZL2,RHO,
     *          IFAULT)
            P = MIN(1.0D0,MAX(P,0.0D0))
         ELSE IF (IND.EQ.1) THEN
            ZU1 = (U1-AM1)/SD1
            ZU2 = (U2-AM2)/SD2
            IFAULT = 0
            P = G01HAF(ZU1,ZU2,RHO,IFAULT)
         ELSE
            ZL1 = (L1-AM1)/SD1
            ZL2 = (L2-AM2)/SD2
            IFAULT = 0
            P = G01HAF(-ZL1,-ZL2,RHO,IFAULT)
         END IF
         G01HBZ = EXP(E)*P
      ELSE
         G01HBZ = 0.0D0
      END IF
      RETURN
      END
