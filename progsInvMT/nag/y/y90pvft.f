      SUBROUTINE Y90PVF(MATRIX,LINE,M,N,A,IA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90PVF :  Print an Integer Matrix  *
C         =======================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IA, M, N
      CHARACTER*1       MATRIX
      CHARACTER*(*)     LINE
C     .. Array Arguments ..
      INTEGER           A(IA,*)
C     .. Local Scalars ..
      INTEGER           I, IREC, J, N1, N2, NP, NREC, UNIT1
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
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
C
      DO 40 I = 1, M
         CALL X04BAF(UNIT1,' ')
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
         NREC = (NP-N2+2)/3
         DO 20 IREC = 1, NREC
            N1 = N2 + 1
            N2 = MIN(N2+3,NP)
            WRITE (REC,FMT=99999) (I,J,A(I,J),J=N1,N2)
            IF (N2-N1.NE.2) REC(26+(N2-N1)*25:) = ' '
            CALL X04BAF(UNIT1,REC)
   20    CONTINUE
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90PVF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (1X,1P,3(5X,'(',I3,',',I3,') ',I8))
      END
