      SUBROUTINE G07BBZ(XMU,XSIG,N,TOLM,MAXITS,T,S,L11,L12,L22,K,WK,
     *                  IERROR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     ROUTINE FOR COMPUTING MAXIMUM-LIKELIHOOD ESTIMATES
C     FOR PARAMETERS OF THE NORMAL DISTRIBUTION FROM
C     GROUPED AND CENSORED DATA USING THE NEWTON-RAPHSON ALGORITHM.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  L11, L12, L22, S, T, TOLM, XMU, XSIG
      INTEGER           IERROR, K, MAXITS
C     .. Array Arguments ..
      DOUBLE PRECISION  WK(*)
      INTEGER           N(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  DS, DXMU, DXSIG, L1, L2, RELMU, RELSIG, TEMP
      INTEGER           I, Q
C     .. External Subroutines ..
      EXTERNAL          G07BBX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      Q = -1
      K = 0
      DS = 0.0D0
C
C     START ITERATIONS
C
      DO 20 I = 1, MAXITS
         K = K + 1
         CALL G07BBX(XMU,XSIG,N,T,S,L1,L2,L11,L12,L22,WK)
         TEMP = (L11*L22-L12*L12)
         IF (TEMP.EQ.0.0D0) THEN
            GO TO 60
         ELSE
            DXSIG = (L1*L12-L2*L11)/TEMP
            DXMU = (-L1*L22+L2*L12)/TEMP
            XMU = XMU + DXMU
            IF (XMU.EQ.0.0D0) THEN
               RELMU = ABS(DXMU)
            ELSE
               RELMU = ABS(DXMU/XMU)
            END IF
            XSIG = XSIG + DXSIG
            IF (ABS(DXSIG).GT.DS) THEN
               Q = Q + 1
            ELSE
               Q = 0
            END IF
            DS = ABS(DXSIG)
            IF (XSIG.LE.0.0D0) XSIG = 0.5D0*(XSIG-DXSIG)
            RELSIG = ABS(DXSIG/XSIG)
            IF (Q.GE.3) THEN
               GO TO 60
            ELSE IF (RELMU.LE.TOLM .AND. RELSIG.LE.TOLM) THEN
               GO TO 40
            END IF
         END IF
   20 CONTINUE
      IERROR = 2
   40 CALL G07BBX(XMU,XSIG,N,T,S,L1,L2,L11,L12,L22,WK)
      RETURN
   60 IERROR = 3
      RETURN
      END
