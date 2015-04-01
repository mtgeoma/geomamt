      SUBROUTINE Y90PDF(MATRIX,LINE,M,N,A,LDA)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ===================================
C         *  Y90PDF :  Print a Real Matrix  *
C         ===================================
C
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           LDA, M, N
      CHARACTER*1       MATRIX
      CHARACTER*(*)     LINE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DUMMY
      INTEGER           I, IREC, J, N1, N2, NP, NREC, UNIT1
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
      CALL X04BAF(UNIT1,' ')
      CALL X04BAF(UNIT1,REC)
C
      IF (X02BEF(DUMMY).LE.10) THEN
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
               REC(28+(N2-N1)*26:) = ' '
               CALL X04BAF(UNIT1,REC)
   20       CONTINUE
   40    CONTINUE
C
      ELSE
         DO 80 I = 1, M
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
            NREC = (NP-N2+1)/2
            DO 60 IREC = 1, NREC
               N1 = N2 + 1
               N2 = MIN(N2+2,NP)
               WRITE (REC,FMT=99998) (I,J,A(I,J),J=N1,N2)
               IF (N1.EQ.N2) REC(38:) = ' '
               CALL X04BAF(UNIT1,REC)
   60       CONTINUE
   80    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90PDF
C
C-----------------------------------------------------------------------
      RETURN
C
99999 FORMAT (1P,3(3X,'(',I3,',',I3,') ',D13.6))
99998 FORMAT (4X,1P,2('(',I3,',',I3,')  ',D19.12,8X))
      END
