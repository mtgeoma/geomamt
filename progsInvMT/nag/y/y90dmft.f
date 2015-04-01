      SUBROUTINE Y90DMF(VMAN,VEXP,SCALE)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==============================
C         *  Y90DMF :  Binary Scaling  *
C         ==============================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      COMPLEX*16        CZERO, CONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,
     *                  CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
C     .. Scalar Arguments ..
      COMPLEX*16        VMAN
      INTEGER           SCALE, VEXP
C     .. Local Scalars ..
      COMPLEX*16        LOWB, RSCALE
      LOGICAL           LOOP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCMPLX, DBLE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Processing
C
C-----------------------------------------------------------------------
      RSCALE = DCMPLX(TWO**SCALE,ZERO)
      LOWB = CONE/RSCALE
C
      IF (VMAN.NE.CZERO) THEN
   20    CONTINUE
         IF (ABS(VMAN).GE.ONE) THEN
            VMAN = VMAN*LOWB
            VEXP = VEXP + SCALE
            LOOP = .TRUE.
         ELSE
            LOOP = .FALSE.
         END IF
         IF (LOOP) GO TO 20
C
   40    CONTINUE
         IF (ABS(VMAN).LT.DBLE(LOWB)) THEN
            VMAN = VMAN*RSCALE
            VEXP = VEXP - SCALE
            LOOP = .TRUE.
         ELSE
            LOOP = .FALSE.
         END IF
         IF (LOOP) GO TO 40
C
      ELSE
         VMAN = CZERO
         VEXP = 0
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90DMF
C
C-----------------------------------------------------------------------
      RETURN
      END
