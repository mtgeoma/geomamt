      SUBROUTINE Y90PXF(LINE,NVECCH,VECCH,IVECCH)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ========================================
C         *  Y90PXF :  Print a Character Vector  *
C         ========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECCH, NVECCH
      CHARACTER*(*)     LINE
C     .. Array Arguments ..
      CHARACTER*(*)     VECCH(*)
C     .. Local Scalars ..
      INTEGER           I, IREC, J, K, L, N1, N2, NREC, UNIT1
      CHARACTER*80      REC
C     .. Local Arrays ..
      CHARACTER*30      ITEM(2)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MIN, MOD
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Print
C
C-----------------------------------------------------------------------
      CALL X04AAF(0,UNIT1)
C
      REC = ' '
      REC(5:) = LINE(1:)
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
C
      L = LEN(VECCH(1))
C
      NREC = (NVECCH+1)/2
      N2 = 0
      DO 40 IREC = 1, NREC
         N1 = N2 + 1
         N2 = MIN(N2+2,NVECCH)
         DO 20 I = N1, N2
            K = MOD(I+1,2) + 1
            J = (I-1)*IVECCH + 1
            ITEM(K) = '"'//VECCH(J)//'"'
            IF (L.GT.28) ITEM(K) (27:) = '..."'
   20    CONTINUE
         WRITE (REC,FMT=99999) (I,ITEM(MOD(I+1,2)+1),I=N1,N2)
         IF (N2-N1.LE.0) REC(43:) = ' '
         CALL X04BAF(UNIT1,REC)
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PXF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,2('(',I3,')',1X,A30,3X))
      END
