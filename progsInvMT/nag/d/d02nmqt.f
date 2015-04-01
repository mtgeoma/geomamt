      SUBROUTINE D02NMQ(WM,IWM,X,N,IER)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 15 REVISED. IER-900 (APR 1991). (LAPACK)
C
C     OLD NAME SOLSF
C
C     .. Scalar Arguments ..
      INTEGER           IER, N
C     .. Array Arguments ..
      DOUBLE PRECISION  WM(*), X(*)
      INTEGER           IWM(*)
C     .. Scalars in Common ..
      INTEGER           IWP
      CHARACTER*6       SINGLR
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, F07AEG
C     .. Common blocks ..
      COMMON            /MD02NM/SINGLR
      COMMON            /SD02NM/IWP
C     .. Save statement ..
      SAVE              /MD02NM/, /SD02NM/
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C THIS ROUTINE MANAGES THE SOLUTION OF THE FULL LINEAR SYSTEM ARISING
C FROM A CHORD ITERATION .
C COMMUNICATION WITH SOLSY USES THE FOLLOWING VARIABLES..
C WM    = REAL WORK SPACE CONTAINING LU DECOMPOSITION OF THE MATRIX
C         STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
C IWM   = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING AT
C         IWM(1).
C X     = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION VECTOR
C         ON OUTPUT, OF LENGTH N.
C IER   = OUTPUT FLAG .  IER = 0 IF NO TROUBLE OCCURRED.
C                        IER = -1 OTHERWISE.
C-----------------------------------------------------------------------
      IF (SINGLR.NE.'NSING2') THEN
         CALL D02NNQ(
     *' BACK SUBSTITUTION TRIED
     *  WITH EITHER SINGULAR MATRIX OR WHEN THE JACOBIAN
     *  HAS NOT BEEN FORMED OR THE ARRAY RWORK HAS BEEN
     *  OVERWRITTEN',1,0,0,0,0,0.0D0,0.0D0)
         IER = -1
         RETURN
      END IF
      CALL F07AEG('No transpose',N,1,WM,N,WM(IWP),X,N,IER)
      RETURN
      END
