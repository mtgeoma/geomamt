      SUBROUTINE G13DJV(I,LMAX,IK,TRAN,PREDZ,SEFZ,LLL,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           I, IFAULT, IK, LLL, LMAX
      CHARACTER*1       TRAN
C     .. Array Arguments ..
      DOUBLE PRECISION  PREDZ(IK,LMAX), SEFZ(IK,LMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION  LOWLIM, TM, TP, VL, WL
      INTEGER           L
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, SQRT
C     .. Executable Statements ..
C
C     correct predictions and their standard errors
C
      IFAULT = 0
      LOWLIM = LOG(X02AMF())
      IF (TRAN.EQ.'L' .OR. TRAN.EQ.'l') THEN
         DO 20 L = LLL, LMAX
            WL = PREDZ(I,L) + 0.5D0*SEFZ(I,L)
            TM = SEFZ(I,L)
            IF (TM.GT.-LOWLIM) THEN
               IFAULT = 1
               RETURN
            ELSE
               TM = EXP(TM) - 1.0D0
            END IF
            VL = 2.0D0*WL
            IF (VL.GT.-LOWLIM) THEN
               IFAULT = 1
               RETURN
            ELSE IF (VL.LT.LOWLIM) THEN
               VL = 0.0D0
            ELSE
               VL = EXP(VL)
            END IF
C
            IF (TM.GE.0.0D0) THEN
               SEFZ(I,L) = SQRT(VL)*SQRT(TM)
            ELSE
               IFAULT = 2
               RETURN
            END IF
C
            IF (WL.GT.-LOWLIM) THEN
               IFAULT = 1
               RETURN
            ELSE IF (WL.LT.LOWLIM) THEN
               PREDZ(I,L) = 0.0D0
            ELSE
               PREDZ(I,L) = EXP(WL)
            END IF
   20    CONTINUE
      ELSE IF (TRAN.EQ.'S' .OR. TRAN.EQ.'s') THEN
         DO 40 L = LLL, LMAX
            TP = PREDZ(I,L)
            TM = SEFZ(I,L)
            WL = TP*TP + TM
            VL = 4.0D0*TP*TP*TM + 2.0D0*TM*TM
            PREDZ(I,L) = WL
            IF (VL.GE.0.0D0) THEN
               SEFZ(I,L) = SQRT(VL)
            ELSE
               IFAULT = 2
               RETURN
            END IF
   40    CONTINUE
      ELSE
         DO 60 L = LLL, LMAX
            IF (SEFZ(I,L).GE.0.0D0) THEN
               SEFZ(I,L) = SQRT(SEFZ(I,L))
            ELSE
               IFAULT = 2
               RETURN
            END IF
   60    CONTINUE
      END IF
C
C
      RETURN
      END
