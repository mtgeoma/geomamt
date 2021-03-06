      SUBROUTINE D02NCF(NEQ,NEQMAX,T,TOUT,Y,YDOTI,RWORK,RTOL,ATOL,ITOL,
     *                  INFORM,FCN,YSAVE,NY2DIM,JAC,WKJAC,NWKJAC,JACPVT,
     *                  NJCPVT,MONITR,ITASK,ITRACE,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NCF')
      INTEGER           NQU, HU, H, HMIN, HMAX, TN, EL0
      PARAMETER         (NQU=10,HU=15,H=16,HMIN=17,HMAX=18,TN=19,EL0=20)
      INTEGER           ML, MU
      PARAMETER         (ML=37,MU=36)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TOUT
      INTEGER           IFAIL, ITASK, ITOL, ITRACE, NEQ, NEQMAX, NJCPVT,
     *                  NWKJAC, NY2DIM
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), RTOL(*), RWORK(4*NEQMAX+50),
     *                  WKJAC(NWKJAC), Y(NEQMAX), YDOTI(NEQMAX),
     *                  YSAVE(NEQMAX,NY2DIM)
      INTEGER           INFORM(23), JACPVT(NJCPVT)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, JAC, MONITR
C     .. Scalars in Common ..
      CHARACTER*6       LAOPTN
C     .. Local Scalars ..
      INTEGER           IFAIL1, IMON, INLN, IRES, IREVCM, LACOR, LEWT,
     *                  LSAVR
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NMF, D02NNQ
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Common blocks ..
      COMMON            /AD02NC/LAOPTN
C     .. Save statement ..
      SAVE              /AD02NC/
C     .. Executable Statements ..
      IF (LAOPTN.NE.'BANDED') THEN
         IFAIL1 = 15
         GO TO 180
      END IF
C
      LEWT = 51
      LACOR = LEWT + NEQMAX
      LSAVR = LACOR + NEQMAX
      IREVCM = 0
      IFAIL1 = 1
   20 CALL D02NMF(NEQ,NEQMAX,T,TOUT,Y,YDOTI,RWORK,RTOL,ATOL,ITOL,INFORM,
     *            YSAVE,NY2DIM,WKJAC,NWKJAC,JACPVT,NJCPVT,IMON,INLN,
     *            IRES,IREVCM,ITASK,ITRACE,IFAIL1)
C
      GO TO (40,100,40,60,80,100,100,120,140,
     *       160) IREVCM
      GO TO 180
   40 CALL FCN(NEQ,RWORK(TN),Y,RWORK(LSAVR),IRES)
      GO TO 20
   60 CALL FCN(NEQ,RWORK(TN),Y,RWORK(LACOR),IRES)
      GO TO 20
   80 CALL FCN(NEQ,RWORK(TN),Y,YDOTI,IRES)
      GO TO 20
  100 CALL D02NNQ(
     *' IMPOSSIBLE VALUE OF THE REVERSE                          COMMUNI
     *CATION VARIABLE ',1,0,0,0,0,0.0D0,0.0D0)
      GO TO 180
  120 CALL JAC(NEQ,RWORK(TN),Y,RWORK(H),RWORK(EL0),INT(RWORK(ML)),
     *         INT(RWORK(MU)),WKJAC)
      GO TO 20
  140 CALL MONITR(NEQ,NEQMAX,RWORK(TN),RWORK(HU),RWORK(H),Y,YDOTI,YSAVE,
     *            RWORK(LSAVR),RWORK(LEWT),IMON,INLN,RWORK(HMIN),
     *            RWORK(HMAX),INT(RWORK(NQU)))
      GO TO 20
  160 CONTINUE
      GO TO 20
  180 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
