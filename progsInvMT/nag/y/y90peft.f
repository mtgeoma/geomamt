      SUBROUTINE Y90PEF(MATRIX,LINE,MNAME1,MNAME2,M,N,A,IA,B,IB)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =====================================
C         *  Y90PEF :  Print Two Real Matrix  *
C         =====================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IA, IB, M, N
      CHARACTER*1       MATRIX
      CHARACTER*(*)     LINE, MNAME1, MNAME2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), B(IB,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DUMMY
      INTEGER           I, J, J1, J2, UNIT1
      CHARACTER*10      NAME1, NAME2
      CHARACTER*80      REC
C     .. External Functions ..
      INTEGER           X02BEF
      LOGICAL           Y90WAF
      EXTERNAL          X02BEF, Y90WAF
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
      NAME1(1:) = MNAME1
      NAME2(1:) = MNAME2
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
C
      IF (X02BEF(DUMMY).LE.10) THEN
         WRITE (REC,FMT=99999) NAME1, NAME2
      ELSE
         WRITE (REC,FMT=99997) NAME1, NAME2
      END IF
      CALL X04BAF(UNIT1,REC)
C
      DO 60 I = 1, M
         IF (Y90WAF(MATRIX,'L') .OR. Y90WAF(MATRIX,'Y')) THEN
            J1 = 1
            J2 = MIN(I,N)
         ELSE IF (Y90WAF(MATRIX,'U') .OR. Y90WAF(MATRIX,'S')) THEN
            J1 = I
            J2 = N
         ELSE
            J1 = 1
            J2 = N
         END IF
         CALL X04BAF(UNIT1,' ')
C
         IF (X02BEF(DUMMY).LE.10) THEN
            DO 20 J = J1, J2
               WRITE (REC,FMT=99998) I, J, A(I,J), B(I,J)
               CALL X04BAF(UNIT1,REC)
   20       CONTINUE
         ELSE
            DO 40 J = J1, J2
               WRITE (REC,FMT=99996) I, J, A(I,J), B(I,J)
               CALL X04BAF(UNIT1,REC)
   40       CONTINUE
         END IF
   60 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PEF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (21X,A10,13X,A10)
99998 FORMAT (5X,'(',I3,',',I3,')',1P,6X,D13.6,10X,D13.6)
99997 FORMAT (24X,A10,19X,A10)
99996 FORMAT (5X,'(',I3,',',I3,')',1P,6X,D19.12,10X,D19.12)
      END
