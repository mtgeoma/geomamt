      SUBROUTINE G04BCY(IBLOCK,KBLOCK,NT,IT,C,LDC,CONST)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C      Computes the NN' matrix
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONST
      INTEGER           IBLOCK, KBLOCK, LDC, NT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,NT)
      INTEGER           IT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  SCALE
      INTEGER           I, II, INCI, INCJ, J, JJ, JTH, K, KK, KTH,
     *                  NBLOCK
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      NBLOCK = ABS(IBLOCK)
      IF (IBLOCK.LT.0) THEN
         INCI = 1
         INCJ = NBLOCK
      ELSE
         INCI = KBLOCK
         INCJ = 1
      END IF
      SCALE = CONST/DBLE(KBLOCK)
      II = 1
      DO 80 I = 1, NBLOCK
         JJ = II
         DO 60 J = 1, KBLOCK
            JTH = IT(JJ)
            KK = II
            DO 20 K = 1, J - 1
               IF (JTH.EQ.IT(KK)) THEN
                  C(JTH,JTH) = C(JTH,JTH) - SCALE
               END IF
               KK = KK + INCJ
   20       CONTINUE
            C(JTH,JTH) = C(JTH,JTH) - SCALE
            KK = KK + INCJ
            DO 40 K = J + 1, KBLOCK
               KTH = IT(KK)
               IF (JTH.LE.KTH) THEN
                  C(KTH,JTH) = C(KTH,JTH) - SCALE
               ELSE
                  C(JTH,KTH) = C(JTH,KTH) - SCALE
               END IF
               KK = KK + INCJ
   40       CONTINUE
            JJ = JJ + INCJ
   60    CONTINUE
         II = II + INCI
   80 CONTINUE
      RETURN
      END
