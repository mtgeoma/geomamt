      SUBROUTINE Y90RGF(DET,VTYPE,N,V,VBOUND,COND,SCALE,DETMAN,DETEXP,
     *                  DIST,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ======================================
C         *  Y90RGF :  Random Vector Generator  *
C         ======================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DETMAN, SCALE
      INTEGER           DETEXP, DIST, N, VTYPE
      CHARACTER*1       DET
C     .. Array Arguments ..
      DOUBLE PRECISION  V(*), VBOUND(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I
C     .. External Functions ..
      DOUBLE PRECISION  Y90TBF
      LOGICAL           Y90WAF
      EXTERNAL          Y90TBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          Y90SMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Generate vector
C
C-----------------------------------------------------------------------
      IF (ABS(VTYPE).EQ.1) THEN
         DO 20 I = 1, N
            IF (DIST.LE.1) THEN
               V(I) = VBOUND(1) + (VBOUND(2)-VBOUND(1))*Y90TBF(DIST,
     *                SEED)
            ELSE IF (DIST.EQ.2) THEN
               V(I) = (VBOUND(1)+VBOUND(2)+(VBOUND(2)-VBOUND(1))
     *                *Y90TBF(DIST,SEED))*HALF
            ELSE
               V(I) = VBOUND(2)*Y90TBF(DIST,SEED) + VBOUND(1)
            END IF
   20    CONTINUE
      ELSE IF (ABS(VTYPE).EQ.2) THEN
         V(1) = ONE
         DO 40 I = 2, N
            V(I) = COND**(-DBLE(I-1)/DBLE(N-1))
   40    CONTINUE
      ELSE IF (ABS(VTYPE).EQ.3) THEN
         V(1) = ONE
         DO 60 I = 2, N
            V(I) = DBLE(1-DBLE(I-1)/DBLE(N-1)) + (DBLE(I-1)/DBLE(N-1))
     *             /COND
   60    CONTINUE
      END IF
C
      IF (ABS(VTYPE).NE.1) THEN
         DO 80 I = 1, N
            V(I) = V(I)*SCALE
   80    CONTINUE
      END IF
C
      IF (VTYPE.LE.-2) THEN
         DO 100 I = 1, N/2
            TEMP = V(I)
            V(I) = V(N-I+1)
            V(N-I+1) = TEMP
  100    CONTINUE
      END IF
C
C     Calculate the determinant
C
      IF (Y90WAF(DET,'D')) THEN
         DETMAN = ONE
         DETEXP = 0
         DO 120 I = 1, N
            DETMAN = DETMAN*V(I)
            CALL Y90SMF(DETMAN,DETEXP,4)
  120    CONTINUE
      ELSE
         DETMAN = ZERO
         DETEXP = 0
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90RGF
C
C-----------------------------------------------------------------------
      RETURN
      END
