      SUBROUTINE D03EDT(U,LEV,NGP,NGRIDX,NGRIDY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Output of vector U on GRID(LEV)
C     (at most NN numbers per line)
C
C     .. Parameters ..
      INTEGER           NN
      PARAMETER         (NN=10)
C     .. Scalar Arguments ..
      INTEGER           LEV
C     .. Array Arguments ..
      DOUBLE PRECISION  U(*)
      INTEGER           NGP(12), NGRIDX(12), NGRIDY(12)
C     .. Local Scalars ..
      INTEGER           I, J, JEND, K, KEND, L, M, N, NADV, NP
      CHARACTER*120     REC
C     .. External Subroutines ..
      EXTERNAL          X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
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
         WRITE (REC,FMT=99998) J
         CALL X04BAF(NADV,REC)
         NP = NP + NGRIDX(LEV)
         L = (KEND-1)/NN
         DO 20 I = 1, L + 1
            M = (I-1)*NN + 1
            N = MIN(KEND,M+NN-1)
            WRITE (REC,FMT=99997) M, N, (U(NP+K),K=M,N)
            CALL X04BAF(NADV,REC)
   20    CONTINUE
   40 CONTINUE
      RETURN
C
99999 FORMAT (' Multigrid Level = ',I6)
99998 FORMAT (' Y-index=',I3)
99997 FORMAT ('  X-index=',I3,'-',I3,'  ',1P,10D10.2)
      END
