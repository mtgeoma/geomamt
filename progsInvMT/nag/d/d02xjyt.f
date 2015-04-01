      SUBROUTINE D02XJY(T,K,YH,NYH,DKY,IFLAG,NEQ,H,TN,HU,NQ,ODCODE)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C-----------------------------------------------------------------------
C  INTERPOLATION ROUTINE DRIVER FOR BLEND GEAR THETA AND DASSL IN D02N
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HU, ODCODE, T, TN
      INTEGER           IFLAG, K, NEQ, NQ, NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  DKY(*), YH(NYH,*)
C     .. External Subroutines ..
      EXTERNAL          D02MVW, D02MWU, D02XJZ
C     .. Executable Statements ..
      IF (ODCODE.LE.2.0D0) THEN
C           SPGEAR OR SBLEND INTERPOLANT
         CALL D02XJZ(T,K,YH,NYH,DKY,IFLAG,NEQ,H,TN,HU,NQ)
      ELSE IF (ODCODE.EQ.3.0D0) THEN
C           STHETB INTERPOLANT
         CALL D02MWU(T,K,YH,NYH,DKY,IFLAG,NEQ,H,TN,HU)
      ELSE IF (ODCODE.EQ.4.0D0) THEN
C           SPDASL INTERPOLANT
         CALL D02MVW(TN,T,DKY,NEQ,NQ,YH,NYH,K)
      END IF
      RETURN
C----------------------- END OF SUBROUTINE D02XJY ----------------------
      END
