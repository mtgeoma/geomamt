      SUBROUTINE G13DPX(NY,ICOL,N,IP,Q,LDQ,RSS,B,SE,COV,TOL,WK,ILLCON)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     COMPUTES REGRESSION PARAMETERS ETC FROM R MATRIX AND Q'Y VECTOR
C     THESE ARE STORED IN Q
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS, TOL
      INTEGER           ICOL, IP, LDQ, N, NY
      LOGICAL           ILLCON
C     .. Array Arguments ..
      DOUBLE PRECISION  B(IP), COV(IP*(IP+1)/2), Q(LDQ,IP+NY), SE(IP),
     *                  WK(IP*(IP+1)/2)
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, RMS
      INTEGER           I, IJ, J
C     .. External Functions ..
      DOUBLE PRECISION  F02WDZ
      EXTERNAL          F02WDZ
C     .. External Subroutines ..
      EXTERNAL          F06EFF, F06PJF, G02AAX, G02AAY, G02AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      IF (ILLCON) THEN
         COND = F02WDZ(IP,Q(1,NY+1),LDQ,WK)
         IF (COND*TOL.LE.1.0D0) ILLCON = .FALSE.
      END IF
      IF ( .NOT. ILLCON) THEN
         CALL F06EFF(IP,Q(1,ICOL),1,B,1)
         CALL F06PJF('U','N','N',IP,Q(1,NY+1),LDQ,B,1)
         IJ = 1
         DO 20 I = 1, IP
            CALL F06EFF(I,Q(1,NY+I),1,COV(IJ),1)
            IJ = IJ + I
   20    CONTINUE
         CALL G02AAZ('U','N',IP,COV)
         CALL G02AAX('U',IP,COV,WK)
         CALL G02AAY('L','N',IP,WK)
         CALL G02AAX('L',IP,WK,COV)
         RMS = RSS/DBLE(N-IP)
         DO 40 I = 1, (IP*IP+IP)/2
            COV(I) = COV(I)*RMS
   40    CONTINUE
         DO 60 J = 1, IP
            SE(J) = SQRT(COV((J*J+J)/2))
   60    CONTINUE
      END IF
C
      RETURN
      END
