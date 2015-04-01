      SUBROUTINE D03EDU(A,LEV,LDA,NGP,NGRIDX,NGRIDY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Output of matrix A on GRID(LEV)
C
C     .. Scalar Arguments ..
      INTEGER           LDA, LEV
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,7)
      INTEGER           NGP(12), NGRIDX(12), NGRIDY(12)
C     .. Local Scalars ..
      INTEGER           I, J, JEND, K, KEND, NADV, NP
      CHARACTER*120     REC
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF
C     .. Executable Statements ..
      CALL X04ABF(0,NADV)
C
      WRITE (REC,FMT=99999) LEV
      CALL X04BAF(NADV,' ')
      CALL X04BAF(NADV,REC)
      CALL X04BAF(NADV,' ')
      NP = -NGRIDX(LEV) + NGP(LEV) - NGRIDX(LEV)*NGRIDY(LEV)
      JEND = NGRIDY(LEV)
      KEND = NGRIDX(LEV)
      DO 40 J = 1, JEND
         NP = NP + NGRIDX(LEV)
         WRITE (REC,FMT=99998) J
         CALL X04BAF(NADV,REC)
         DO 20 K = 1, KEND
            WRITE (REC,FMT=99997) K, (A(NP+K,I),I=1,7)
            CALL X04BAF(NADV,REC)
   20    CONTINUE
   40 CONTINUE
      RETURN
C
99999 FORMAT (' Multigrid Level =',I6)
99998 FORMAT (' Y-index=',I3)
99997 FORMAT ('  X-index=',I3,'  ',1P,7D10.2)
      END
