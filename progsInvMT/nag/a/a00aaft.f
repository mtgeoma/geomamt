      SUBROUTINE A00AAF
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Writes information about the particular implementation of the
C     NAG Library in use.
C
C     The output channel is given by a call to X04ABF.
C
C     .. Local Scalars ..
      INTEGER          I, NADV
C     .. Local Arrays ..
      CHARACTER*80     MSG(20)
C     .. External Subroutines ..
      EXTERNAL         A00AAZ, X04ABF, X04BAF
C     .. Executable Statements ..
      CALL A00AAZ(MSG)
      CALL X04ABF(0,NADV)
      DO 20 I = 1, 7
         CALL X04BAF(NADV,MSG(I))
   20 CONTINUE
      CALL X04BAF(NADV,MSG(20))
      RETURN
      END
