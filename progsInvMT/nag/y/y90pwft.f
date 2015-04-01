      SUBROUTINE Y90PWF(MATRIX,LINE,MNAME1,MNAME2,M,N,A,IA,B,IB)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==========================================
C         *  Y90PWF :  Print Two Integer Matrices  *
C         ==========================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IA, IB, M, N
      CHARACTER*1       MATRIX
      CHARACTER*(*)     LINE, MNAME1, MNAME2
C     .. Array Arguments ..
      INTEGER           A(IA,*), B(IB,*)
C     .. Local Scalars ..
      INTEGER           I, J, K, N1, N2, NP, NREC, UNIT1
      CHARACTER*10      NAME1, NAME2
      CHARACTER*80      REC
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
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
C
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
      IF (N.LE.1) THEN
         WRITE (REC,FMT=99999) NAME1, NAME2
      ELSE
         WRITE (REC,FMT=99999) NAME1, NAME2, NAME1, NAME2
      END IF
      CALL X04BAF(UNIT1,REC)
C
      DO 40 I = 1, M
         IF (Y90WAF(MATRIX,'L') .OR. Y90WAF(MATRIX,'Y')) THEN
            N2 = 0
            NP = MIN(I,N)
         ELSE IF (Y90WAF(MATRIX,'U') .OR. Y90WAF(MATRIX,'S')) THEN
            N2 = I - 1
            NP = N
         ELSE
            N2 = 0
            NP = N
         END IF
         CALL X04BAF(UNIT1,' ')
         NREC = (NP-N2+1)/2
         DO 20 J = 1, NREC
            N1 = N2 + 1
            N2 = MIN(N2+2,NP)
            WRITE (REC,FMT=99998) (I,K,A(I,K),B(I,K),K=N1,N2)
            IF (N2.LE.N1) REC(38:) = ' '
            CALL X04BAF(UNIT1,REC)
   20    CONTINUE
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PWF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (5X,2(9X,A10,2X,A10,5X))
99998 FORMAT (5X,2('(',I3,',',I3,')',1P,2X,I8,4X,I8,5X))
      END
