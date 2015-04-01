      SUBROUTINE Y90PCF(MATRIX,LINE,MNAME1,MNAME2,M,N,A,IA,B,IB)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==========================================
C         *  Y90PCF :  Print Two Complex Matrices  *
C         ==========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IA, IB, M, N
      CHARACTER*1       MATRIX
      CHARACTER*(*)     LINE, MNAME1, MNAME2
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), B(IB,*)
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
      INTRINSIC         DIMAG, MIN, DBLE
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
C
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
      IF (X02BEF(DUMMY).LE.10) THEN
         WRITE (REC,FMT=99999) NAME1, NAME2
      ELSE
         WRITE (REC,FMT=99997) NAME1, NAME2
      END IF
      CALL X04BAF(UNIT1,REC)
C
      DO 60 I = 1, M
         IF (Y90WAF(MATRIX,'L') .OR. Y90WAF(MATRIX,'E')) THEN
            J1 = 1
            J2 = MIN(I,N)
         ELSE IF (Y90WAF(MATRIX,'U') .OR. Y90WAF(MATRIX,'H')) THEN
            J1 = I
            J2 = N
         ELSE
            J1 = 1
            J2 = N
         END IF
         CALL X04BAF(UNIT1,' ')
         IF (X02BEF(DUMMY).LE.10) THEN
            DO 20 J = J1, J2
               WRITE (REC,FMT=99998) I, J, A(I,J), B(I,J)
               CALL X04BAF(UNIT1,REC)
   20       CONTINUE
         ELSE
            DO 40 J = J1, J2
               WRITE (REC,FMT=99996) I, J, DBLE(A(I,J)), DBLE(B(I,J))
               CALL X04BAF(UNIT1,REC)
               WRITE (REC,FMT=99995) DIMAG(A(I,J)), DIMAG(B(I,J))
               CALL X04BAF(UNIT1,REC)
   40       CONTINUE
         END IF
   60 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PCF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (25X,A10,25X,A10)
99998 FORMAT (5X,'(',I3,',',I3,')',1P,3X,D13.6,2X,D13.6,7X,D13.6,2X,
     *       D13.6)
99997 FORMAT (24X,A10,19X,A10)
99996 FORMAT (4X,'(',I3,',',I3,')',1P,6X,D19.12,10X,D19.12,'  REAL')
99995 FORMAT (13X,1P,6X,D19.12,10X,D19.12,'  IMAG')
      END
