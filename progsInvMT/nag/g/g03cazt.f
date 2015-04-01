      SUBROUTINE G03CAZ(K,X,F,G,ISTATE,GP,COND,PD,NIT,NF,IW,LIW,W,LW)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     prints iteration information
C     for maximum likelihood Factor Analysis
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, F, GP
      INTEGER           K, LIW, LW, NF, NIT
      LOGICAL           PD
C     .. Array Arguments ..
      DOUBLE PRECISION  G(K), W(LW), X(K)
      INTEGER           ISTATE(K), IW(LIW)
C     .. Local Scalars ..
      DOUBLE PRECISION  COMM
      INTEGER           I, IVAR, NADV
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF
C     .. Executable Statements ..
      CALL X04ABF(0,NADV)
      WRITE (REC,FMT=99999) NIT, NF
      CALL X04BAF(NADV,REC(1))
      CALL X04BAF(NADV,REC(2))
      WRITE (REC,FMT=99993) F
      CALL X04BAF(NADV,REC(1))
      WRITE (REC,FMT=99998)
      CALL X04BAF(NADV,REC(1))
      CALL X04BAF(NADV,REC(2))
      CALL X04BAF(NADV,REC(3))
      IVAR = IW(4) + 2*K*K + K
      DO 20 I = 1, K
         COMM = (W(IVAR+I)-X(I))/W(IVAR+I)
         IF (ISTATE(I).GE.0) THEN
            WRITE (REC,FMT=99997) I, COMM
         ELSE IF (ISTATE(I).EQ.-1) THEN
            WRITE (REC,FMT=99996) I, COMM
         ELSE IF (ISTATE(I).EQ.-2) THEN
            WRITE (REC,FMT=99995) I, COMM
         ELSE
            WRITE (REC,FMT=99994) I, COMM
         END IF
         CALL X04BAF(NADV,REC(1))
   20 CONTINUE
      RETURN
C
99999 FORMAT (/' Iterations performed =',I16,', function evaluations =',
     *       I16)
99998 FORMAT (/13X,'Variable',4X,'Standardized',/25X,'Communalities')
99997 FORMAT (1X,I16,11X,F7.4)
99996 FORMAT (1X,I16,11X,F7.4,' * at lower bound')
99995 FORMAT (1X,I16,11X,F7.4,' * at upper bound')
99994 FORMAT (1X,I16,11X,F7.4,' * held constant')
99993 FORMAT (' Criterion = ',D16.6)
      END
