      SUBROUTINE G02GBS(IP,IRANK,SVD,Q,LDQ,COV,WK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      INTEGER           IP, IRANK, LDQ
      LOGICAL           SVD
C     .. Array Arguments ..
      DOUBLE PRECISION  COV((IP*IP+IP)/2), Q(LDQ,IP), WK((IP*IP+IP)/2)
C     .. Local Scalars ..
      INTEGER           I, IJ
C     .. External Subroutines ..
      EXTERNAL          G02AAX, G02AAY, G02AAZ, DCOPY, DGEMV
C     .. Executable Statements ..
      IJ = 1
      IF (SVD) THEN
         DO 20 I = 1, IP
            CALL DCOPY(IRANK,Q(1,I),1,WK,1)
            CALL DGEMV('T',IRANK,I,1.0D0,Q,LDQ,WK,1,0.0D0,COV(IJ),1)
            IJ = IJ + I
   20    CONTINUE
      ELSE
         DO 40 I = 1, IP
            CALL DCOPY(I,Q(1,I),1,COV(IJ),1)
            IJ = IJ + I
   40    CONTINUE
         CALL G02AAZ('U','N',IP,COV)
         CALL G02AAX('U',IP,COV,WK)
         CALL G02AAY('L','N',IP,WK)
         CALL G02AAX('L',IP,WK,COV)
      END IF
      END
