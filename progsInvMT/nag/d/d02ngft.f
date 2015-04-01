      SUBROUTINE D02NGF(NEQ,NEQMAX,T,TOUT,Y,YDOTI,RWORK,RTOL,ATOL,ITOL,
     *                  INFORM,RESID,YSAVE,NY2DIM,JAC,WKJAC,NWKJAC,
     *                  MONITR,LDERIV,ITASK,ITRACE,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      INTEGER           NJCPVT
      PARAMETER         (NJCPVT=1)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02NGF')
      INTEGER           NQU, HU, H, HMIN, HMAX, TN, EL0
      PARAMETER         (NQU=10,HU=15,H=16,HMIN=17,HMAX=18,TN=19,EL0=20)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TOUT
      INTEGER           IFAIL, ITASK, ITOL, ITRACE, NEQ, NEQMAX, NWKJAC,
     *                  NY2DIM
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), RTOL(*), RWORK(4*NEQMAX+50),
     *                  WKJAC(NWKJAC), Y(NEQMAX), YDOTI(NEQMAX),
     *                  YSAVE(NEQMAX,NY2DIM)
      INTEGER           INFORM(23)
      LOGICAL           LDERIV(2)
C     .. Subroutine Arguments ..
      EXTERNAL          JAC, MONITR, RESID
C     .. Scalars in Common ..
      CHARACTER*6       LAOPTN
C     .. Local Scalars ..
      INTEGER           IFAIL1, IMON, INLN, IRES, IREVCM, LACOR, LEWT,
     *                  LSAVR
C     .. Local Arrays ..
      INTEGER           JACPVT(1)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNF
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Common blocks ..
      COMMON            /AD02NC/LAOPTN
C     .. Save statement ..
      SAVE              /AD02NC/
C     .. Executable Statements ..
      IF (LAOPTN.NE.'FULLOP') THEN
         IFAIL1 = 15
         GO TO 180
      END IF
C
      LEWT = 51
      LACOR = LEWT + NEQMAX
      LSAVR = LACOR + NEQMAX
      IREVCM = 0
      IFAIL1 = 1
C
   20 CONTINUE
      CALL D02NNF(NEQ,NEQMAX,T,TOUT,Y,YDOTI,RWORK,RTOL,ATOL,ITOL,INFORM,
     *            YSAVE,NY2DIM,WKJAC,NWKJAC,JACPVT,NJCPVT,IMON,INLN,
     *            IRES,IREVCM,LDERIV,ITASK,ITRACE,IFAIL1)
C
      GO TO (60,40,60,80,100,60,80,120,140,
     *       160,60) IREVCM
C
      GO TO 180
   40 CALL RESID(NEQ,RWORK(TN),Y,RWORK(LACOR),RWORK(LSAVR),IRES)
      GO TO 20
   60 CALL RESID(NEQ,RWORK(TN),Y,YDOTI,RWORK(LSAVR),IRES)
      GO TO 20
   80 CALL RESID(NEQ,RWORK(TN),Y,YDOTI,RWORK(LACOR),IRES)
      GO TO 20
  100 CALL RESID(NEQ,RWORK(TN),Y,RWORK(LSAVR),YDOTI,IRES)
      GO TO 20
  120 CALL JAC(NEQ,RWORK(TN),Y,YDOTI,RWORK(H),RWORK(EL0),WKJAC)
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
