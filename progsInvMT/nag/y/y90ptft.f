      SUBROUTINE Y90PTF(LINE,NVECL,VECL,IVECL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ======================================
C         *  Y90PTF :  Print a Logical Vector  *
C         ======================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IVECL, NVECL
      CHARACTER*(*)     LINE
C     .. Array Arguments ..
      LOGICAL           VECL(*)
C     .. Local Scalars ..
      INTEGER           I, IREC, J, N1, N2, NREC, UNIT1
      CHARACTER*80      REC
C     .. Local Arrays ..
      CHARACTER*8       VALUE(4)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
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
      NREC = (NVECL+3)/4
      N2 = 0
      DO 40 IREC = 1, NREC
         N1 = N2 + 1
         N2 = MIN(N2+4,NVECL)
         DO 20 J = N1, N2
            IF (VECL((J-1)*IVECL+1)) THEN
               VALUE(J-N1+1) = ' .TRUE.'
            ELSE
               VALUE(J-N1+1) = ' .FALSE.'
            END IF
   20    CONTINUE
         WRITE (REC,FMT=99999) (I,VALUE(I-N1+1),I=N1,N2)
         IF (N2-N1.NE.3) REC(24+(N2-N1)*19:) = ' '
         CALL X04BAF(UNIT1,REC)
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PTF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (4X,4('(',I3,')',A8,6X))
      END
