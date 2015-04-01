      SUBROUTINE Y90RRX(PIVOT,NZ,IW1,IW2,IW3,IW4,NK,HIT,LDHIT)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C-----------------------------------------------------------------------
C
C         ===================================================
C         *  Y90RRX :  Utility for Sparse Matrix Generator  *
C         ===================================================
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           LDHIT, NK, NZ, PIVOT
C     .. Array Arguments ..
      INTEGER           IW1(*), IW2(*), IW3(*), IW4(*)
      LOGICAL           HIT(LDHIT,*)
C     .. Local Scalars ..
      INTEGER           I, K
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Find whether any mixing possible at the current level
C
C-----------------------------------------------------------------------
      NK = 0
      I = IW1(NZ)
      IF ((I.NE.PIVOT) .AND. ( .NOT. HIT(I,PIVOT))) THEN
         NK = NK + 1
         IW4(NK) = I
      END IF
      DO 20 K = 1, IW3(NZ) - 1
         I = IW2(I)
         IF ((I.NE.PIVOT) .AND. ( .NOT. HIT(I,PIVOT))) THEN
            NK = NK + 1
            IW4(NK) = I
         END IF
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RRX
C
C-----------------------------------------------------------------------
      RETURN
      END
