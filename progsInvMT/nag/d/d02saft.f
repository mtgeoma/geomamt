      SUBROUTINE D02SAF(P,M,N,N1,PE,PF,E,DP,NPOINT,WP,IWP,ICOUNT,RANGE,
     *                  BC,FCN,EQN,CONSTR,YMAX,MONIT,PRSOL,W,IW1,IW2,
     *                  IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     SOLVES A TWO POINT BVP
C     BC, EQN, FCN, MONIT, PRSOL, RANGE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02SAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  YMAX
      INTEGER           ICOUNT, IFAIL, IW1, IW2, IWP, M, N, N1, NPOINT
C     .. Array Arguments ..
      DOUBLE PRECISION  DP(M), E(N), P(M), PE(M), PF(M), W(IW1,IW2),
     *                  WP(IWP,6)
C     .. Function Arguments ..
      LOGICAL           CONSTR
      EXTERNAL          CONSTR
C     .. Subroutine Arguments ..
      EXTERNAL          BC, EQN, FCN, MONIT, PRSOL, RANGE
C     .. Scalars in Common ..
      DOUBLE PRECISION  DM, EPS, EPSFAC, MACHEP, PNORM, PNORM1, SQEPS
      INTEGER           COUNT, ICASE, IEPS, IFAIL1, IFAIL2, IFAIL3,
     *                  IFLAG1, ISTATE
C     .. Local Scalars ..
      INTEGER           IFLAG, M1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02HAX, D02HAY, D02HAZ, D02HBS, D02SAV
C     .. Common blocks ..
      COMMON            /AD02HB/IFAIL2, IFAIL3
      COMMON            /AD02SA/PNORM, PNORM1, EPS, MACHEP, SQEPS, DM,
     *                  EPSFAC, ISTATE, IFLAG1, COUNT, IFAIL1, IEPS,
     *                  ICASE
C     .. Executable Statements ..
      IFAIL3 = 1
      ICASE = 1
      IFLAG1 = 1
      IFLAG = 1
      M1 = 1
      CALL D02SAV(P,M,W,W,N,N1,PE,PF,E,DP,NPOINT,WP,IWP,WP(1,4),WP(2,4)
     *            ,ICOUNT,D02HAY,RANGE,D02HAZ,BC,D02HBS,D02HAX,FCN,EQN,
     *            CONSTR,YMAX,MONIT,PRSOL,W,M1,W,IW1,IW2,IFLAG)
      IFAIL = P01ABF(IFAIL,IFLAG,SRNAME,0,P01REC)
      RETURN
      END
