      SUBROUTINE D02NMN(WM,IWM,X,N,IER)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-975 (JUN 1993).
C
C     OLD NAME SOLSB
C
C     .. Scalar Arguments ..
      INTEGER           IER, N
C     .. Array Arguments ..
      DOUBLE PRECISION  WM(*), X(*)
      INTEGER           IWM(*)
C     .. Scalars in Common ..
      INTEGER           ICALL, MFILLN, ML, MU
      CHARACTER*6       SINGLR
C     .. Local Scalars ..
      INTEGER           IL, MBAND
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, F04LDF
C     .. Common blocks ..
      COMMON            /MD02NM/SINGLR
      COMMON            /QD02NM/MU, ML, ICALL, MFILLN
C     .. Save statement ..
      SAVE              /QD02NM/, /MD02NM/
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C THIS ROUTINE MANAGES THE SOLUTION OF THE BANDED LINEAR SYSTEM ARISING
C FROM A CHORD ITERATION . IT CALLS THE NAG ROUTINE F04LDF TO PERFORM
C THE BACK SUBSTITUTION .
C COMMUNICATION WITH SOLSY USES THE FOLLOWING VARIABLES..
C WM    = REAL WORK SPACE CONTAINING LU DECOMPOSITION OF THE MATRIX
C IWM   = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION
C X     = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION VECTOR
C         ON OUTPUT, OF LENGTH N.
C IER   = OUTPUT FLAG .  IER = 0 IF NO TROUBLE OCCURRED.
C                        IER = -1 OTHERWISE.
C-----------------------------------------------------------------------
      IER = 1
      IF (SINGLR.NE.'NSING1') THEN
         CALL D02NNQ(
     *' BACK SUBSTITUTION TRIED
     *  WITH SINGULAR JACOBIAN OR BEFORE JACOBIAN HAS BEEN
     *  FORMED OR THE ARRAY RWORK HAS BEEN OVERWRITTEN',1,0,0,0,0,0.0D0,
     *               0.0D0)
         IER = -1
         RETURN
      END IF
      MBAND = ML + MU + 1
      IL = MAX(1,ML)
      CALL F04LDF(N,ML,MU,1,WM,MBAND,WM(MFILLN),IL,IWM,X,N,IER)
      IF (IER.NE.0) IER = -1
      RETURN
      END
