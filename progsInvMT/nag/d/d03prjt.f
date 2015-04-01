      SUBROUTINE D03PRJ(TIME,M,NPDE,NIP,XOP,MAXNPT,CONST,STPRAT,IPMINF,
     *                  FMON,NOP,XNP,DXW,CHISTR,IFDUM)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     MARK 17 REVISED. IER-1554 (JUN 1995).
C ----------------------------------------------------------------------
C
C     Version of remesh for Thornton/Leeds software modules
C     interpolation of solution values onto new mesh is left
C     to calling routine.
C     Remeshing algorithm is that of Kautsky and Nichols
C     using padded monitor function.
C     Calculates a non-uniform mesh where the mesh is adapted
C     to fit the solution behaviour and so is ideally suited
C     for boundary layer type behaviour. The method can be used
C     for 2-pt ODE boundary value problems or for PDE problems
C     with wave/flame/shock fronts.
C
C
C     Adaptive space remeshing routine by R.M.Furzeland (TRC,13/10/83)
C     Method:
C     -------
C    Graded mesh based on padded monitor function of Kautsky and Nichols
C     given the even spaced mesh pts CHISTR(I), graded mesh points
C     XNP(CHISTR) are formed by linear interpolation from the
C     mapping CHI(XOP)=(1/G)*INTEGRAL(M(X)) from X0 to X where the
C     monitor function M(X) is a uniformly distr. variation of F(X)
C     vectors based on monitor function norm (e.g. L2) and weighting
C     and is normalised using G=TOTAL INTEGRAL (M(X))
C     The min. value of M(X) (and hence max delx) is limited by Pereyra
C     concept and, further, the adjacent mesh ratio is controlled by
C     the extra padding of M(X) as described by Kautsky and Nichols
C
C     User parameters:
C     ----------------
C     TIME - current time value (zero if problem does not involve time)
C     NPDE - no. of PDEs or unknowns
C     M - type of geometry
C     NIP -  number of input space mesh points
C   XOP - ARRAY(NIP) of initial space mesh - must be monotonically incr.
C     MAXNPT - maximum 2nd dimension of DU as defined in driver program
C          also used as upper bound for nop (MAXNPT=NPTS in present
C          version)
C     CONST - input bound on integral of monitor over each step. size
C         depends on NIP and monitor. Typical value is 5.0D-1.
C         Controls the local truncation error if the monitor is chosen
C         accordingly
C         Decreasing C increases the mesh ratios, decreasing C then
C         mesh tends to uniform. The user should check that CONST
C         and the max. sub-integral actually produced are approx. equal.
C
C  STPRAT - Input bound on adjacent mesh ratio ( > 1 typically 1.5 to 3)
C          this also controls local truncation errors.
C
C    IPMINF - Controls print out of mesh information. IPMINF=0 no print,
C          =1 summary print of mesh characterisitics e.g. mesh ratios &
C          max. integral of monitor per interval (approx.=input CONST)
C          =2 print of monitors, old & new meshes and mesh sizes
C
C     FMON -  Input array of monitor values evaluated at XOP(I=1,..,NIP)
C         Values are calculated according to the user subroutine MONFFD
C         or MONFKB.
C
C     N.B. D03PRJ scales the input values by integral of FMON .
C    ---- D03PRJ also checks FMON non-negative and not all FMON(I) zero.
C         If IPMINF.GT.0 then on output CHISTR contains monitor
C         values at XNP(I), I=1,...,NOP
C   NOP -   If zero on input then NOP calculated by program else NOP set
C         to NIP by program
C         Output calculated no. of mesh points needed for near optimum
C         set of mesh points XNP
C     XNP - Output array with new mesh values I=1,...,NOP
C     DXW -   ARRAY(NOP) used as work space for monitor G
C         and also for CHI array (see above)
C         On output DXW(I)=mesh intervals XNP(I+1)-XNP(I)  I=1,NOP-1
C     CHISTR(NOP)  - Treat as work space
C         If IPMINF.GT.0 then on output CHISTR contains monitor
C         values at XNP(I), I=1,...,NOP
C     IFDUM - Error flag set to 0 on exit if remesh calculation was ok.
C         Uses nag soft/hard fail principles.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONST, STPRAT, TIME
      INTEGER           IFDUM, IPMINF, M, MAXNPT, NIP, NOP, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  CHISTR(MAXNPT), DXW(MAXNPT), FMON(NIP),
     *                  XNP(MAXNPT), XOP(NIP)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, TOUT, UROUND
      INTEGER           ICOUNT, IDEV, IOVFLO, ITRACE, NRMESH
      LOGICAL           REMESH
C     .. Local Scalars ..
      DOUBLE PRECISION  DRAT, DRMAX, DUN100, FLNP, RINT, RLAM, RMFMAX,
     *                  RMXINT, SUM
      INTEGER           I, I1, IDRMAX, IFAIL, IMODE, IND, ISTAGE, J,
     *                  NIP1, NOM2, NOP1
      CHARACTER*80      REC
      CHARACTER*200     ERRMSG
C     .. External Subroutines ..
      EXTERNAL          D02NNN, D02NNQ, D03PZW, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, LOG, MAX
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /HD03PC/TOUT, REMESH, NRMESH, ICOUNT
C     .. Save statement ..
      SAVE              /AD02NM/, /CD03PC/, /HD03PC/
C     .. Executable Statements ..
C
C     Validates input parameters ..
      ISTAGE = 0
      IFAIL = 0
      DUN100 = 10000.0D0*DUNFLO
C
C     Check validity of input mesh and monitor function.
C     XOP() should be monotonically increasing
C     and FMON should be non-negative and not all zero.
      SUM = 0.0D0
      IF (FMON(1).LT.0.0D0) GO TO 580
      DO 20 I = 2, NIP
         IF (XOP(I).LE.XOP(I-1)) GO TO 580
         IF (FMON(I).LT.0.0D0) GO TO 580
         SUM = SUM + 0.5D0*(FMON(I)+FMON(I-1))*(XOP(I)-XOP(I-1))
   20 CONTINUE
      IF (NRMESH.EQ.0 .AND. TIME.LT.TOUT) THEN
         NOP = NIP
         DO 40 I = 1, NIP
            XNP(I) = XOP(I)
   40    CONTINUE
         IFDUM = 0
         RETURN
      END IF
      IF (SUM.LT.DUN100) THEN
C        Zero monitor function. leave the mesh as it is.
         NOP = NIP
         DO 60 I = 1, NIP
            XNP(I) = XOP(I)
   60    CONTINUE
         IF (ITRACE.GT.0 .OR. IPMINF.GT.0) THEN
            WRITE (REC,FMT=99997)
            CALL X04ABF(0,IDEV)
            CALL X04BAF(IDEV,REC)
         END IF
         IFDUM = 0
         RETURN
      END IF
C
C      Scale input monitor function by integral of monitor function
      DO 80 I = 1, NIP
         FMON(I) = FMON(I)/SUM
   80 CONTINUE
      IF (IPMINF.EQ.0) GO TO 120
      WRITE (REC,FMT=99996) TIME, SUM
      CALL X04ABF(0,IDEV)
      CALL X04BAF(IDEV,REC)
C     Compute max. sub-integral of input scaled monitor
      RMXINT = 0.0D0
      IDRMAX = 1
      DO 100 I = 2, NIP
         RINT = 0.5D0*(FMON(I)+FMON(I-1))*(XOP(I)-XOP(I-1))
         IF (RINT.LE.RMXINT) GO TO 100
         IDRMAX = I
         RMXINT = RINT
  100 CONTINUE
      WRITE (REC,FMT=99981)
      CALL X04BAF(IDEV,REC)
      WRITE (REC,FMT=99979) RMXINT, IDRMAX
      CALL X04BAF(IDEV,REC)
C     Sets mode of working. IMODE=0 is Kautsky and Nichols, IMODE=1 is
C     fixed no. output points NOP=NIP
  120 IMODE = 0
      IF (NOP.EQ.0) GO TO 140
      IMODE = 1
      NOP = NIP
  140 NIP1 = NIP - 1
C-----------------------------------------------------------------------
C     Formation of padded monitor function in the work array DXW ( = G)
C-----------------------------------------------------------------------
C     Forward sweep:
      ISTAGE = 1
      IF (IPMINF.GT.1) THEN
         WRITE (REC,FMT=99998)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(XOP,NIP,29)
         WRITE (REC,FMT=99995)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(FMON,NIP,29)
      END IF
      I = 1
      IND = 0
      RLAM = LOG(STPRAT)/CONST
  160 J = I + 1
  180 IF (J.GT.NIP) GO TO 240
C Padding function DXW(J) fitted at point XOP(I) as a function of XOP(J)
      DXW(J) = FMON(I)/(1.0D0+RLAM*(XOP(J)-XOP(I))*FMON(I))
      IF (DXW(J).LE.FMON(J)) GO TO 200
C     Padding function DXW greater than monitor FMON
      IND = IND + 1
      FMON(J) = DXW(J)
      J = J + 1
      GO TO 180
C     Padding function < monitor function
  200 IF (IND.NE.0) GO TO 220
C     Fit padding function at next point XOP(I) and restart
      I = I + 1
      GO TO 160
C     Completion of padding in interval XOP(I)...XOP(I+IND)
  220 I = J
      IND = 0
      GO TO 160
C     Reverse sweep:
  240 IF (IPMINF.GT.1) THEN
         WRITE (REC,FMT=99994)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(FMON,NIP,29)
      END IF
      I = NIP
      IND = 0
  260 J = I - 1
  280 IF (J.LT.1) GO TO 340
C Padding function DXW(J) fitted at point XOP(I) as a function of XOP(J)
      DXW(J) = FMON(I)/(1.0D0-RLAM*(XOP(J)-XOP(I))*FMON(I))
      IF (DXW(J).LE.FMON(J)) GO TO 300
C     Padding function DXW greater than monitor FMON
      IND = IND + 1
      FMON(J) = DXW(J)
      J = J - 1
      GO TO 280
C     Padding function < monitor function
  300 IF (IND.NE.0) GO TO 320
C     Fit padding function at next point XOP(I) and restart
      I = I - 1
      GO TO 260
C     Completion of padding in interval XOP(I)...XOP(I-IND)
  320 I = J
      IND = 0
      GO TO 260
C-----------------------------------------------------------------------
C Formation of new mesh equi-distributed according to the padded monitor
C     Work array DXW now contains CHI
C-----------------------------------------------------------------------
  340 IF (IPMINF.GT.1) THEN
         WRITE (REC,FMT=99993)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(FMON,NIP,29)
      END IF
      ISTAGE = 2
C     Sum for padded monitor function
      SUM = 0.0D0
      DO 360 I = 2, NIP
         SUM = SUM + 0.5D0*(FMON(I)+FMON(I-1))*(XOP(I)-XOP(I-1))
         DXW(I) = SUM
  360 CONTINUE
      IF (ABS(SUM).GT.DUN100) GO TO 380
      WRITE (REC,FMT=99992)
      CALL X04BAF(IDEV,REC)
      GO TO 580
C     Normalise partial integrals CHI (held in DXW) using SUM
  380 DO 400 I = 2, NIP1
         DXW(I) = DXW(I)/SUM
  400 CONTINUE
      DXW(1) = 0.0D0
      DXW(NIP) = 1.0D0
      IF (IMODE.EQ.1) GO TO 420
C   Calculate variable number of output points for near optimum mesh XNP
C     if IMODE=0 (else if IMODE=1 then fixed NOP=NIP)
      ISTAGE = 3
      NOP = INT(SUM/CONST) + 2
      IF (NOP.GT.MAXNPT) THEN
C     No. output pts. needed exceeds MAXNPT. Sets NOP=MAXNPT & continues
         WRITE (REC,FMT=99999) NOP
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
         NOP = MAXNPT
      END IF
  420 IF (IPMINF.GT.0) THEN
         WRITE (REC,FMT=99991) NIP, NOP
         CALL X04BAF(IDEV,REC)
         WRITE (REC,FMT=99990) SUM, CONST
         CALL X04BAF(IDEV,REC)
      END IF
C ----------------------------------------------------------------------
C     Given transformation CHI(XOP) in DXW(I) then produce graded mesh
C     Corresponding to equi-spaced CHISTR values  (inverse interp.)
C    Linear interp. (not spline) to ensure monoticity of XNP + tests XOP
C ----------------------------------------------------------------------
      ISTAGE = 4
      FLNP = NOP - 1
      DO 440 I = 1, NOP
         CHISTR(I) = (I-1)/FLNP
  440 CONTINUE
C XNP(1)=XOP(1) AND XNP)NOP)=XOP(NIP) ensured by CHISTR end values 0 & 1
      I1 = 1
      CALL D03PZW(CHISTR,XNP,NOP,DXW,M,XOP,NIP,I1,I1,IFAIL)
      IF (IFAIL.NE.0) GO TO 580
      NOP1 = NOP - 1
      DO 460 I = 1, NOP1
         DXW(I) = XNP(I+1) - XNP(I)
         IF (DXW(I).GT.DUN100) GO TO 460
         ERRMSG =
     *' Remeshing has produced zero or negative mesh spacing
     * at point (=I1).  Check size/scaling of monitor function.'
         CALL D02NNQ(ERRMSG,1,1,I,0,0,0.0D0,0.0D0)
         GO TO 580
  460 CONTINUE
      DXW(NOP) = 0.0D0
      IF (IPMINF.EQ.0) GO TO 560
C                              End of main calculation of new mesh
C   --------------------------------------------------------------------
C
C     Print and calculate extra mesh info. if IPMINF>0
C     Max. adjacent mesh ratio
      DRMAX = 0.0D0
      NOM2 = NOP1 - 1
      IDRMAX = 1
      DO 480 I = 1, NOM2
         DRAT = ABS(DXW(I)/DXW(I+1))
         DRAT = MAX(DRAT,1.0D0/DRAT)
         IF (DRAT.LE.DRMAX) GO TO 480
         IDRMAX = I
         DRMAX = DRAT
  480 CONTINUE
C     Ratio of max. mesh to even mesh
      RMFMAX = 0.0D0
      DO 500 I = 1, NOP1
         RMFMAX = MAX(RMFMAX,ABS(DXW(I)))
  500 CONTINUE
      RMFMAX = RMFMAX*FLNP/(XNP(NOP)-XNP(1))
      WRITE (REC,FMT=99989) IDRMAX, DRMAX
      CALL X04ABF(0,IDEV)
      CALL X04BAF(IDEV,REC)
      WRITE (REC,FMT=99988) RMFMAX
      CALL X04BAF(IDEV,REC)
      IF (IPMINF.LT.2) GO TO 520
      WRITE (REC,FMT=99987)
      CALL X04BAF(IDEV,REC)
      CALL D02NNN(XOP,NIP,29)
      WRITE (REC,FMT=99986) NOP
      CALL X04BAF(IDEV,REC)
      CALL D02NNN(XNP,NOP,29)
      WRITE (REC,FMT=99985)
      CALL X04BAF(IDEV,REC)
      CALL D02NNN(DXW,NOP1,29)
C ----------------------------------------------------------------------
C     Calculate monitor on new mesh points and max integral of monitor
C     over each new mesh interval
C ----------------------------------------------------------------------
  520 ISTAGE = 5
C     Linear interpolation to find monitor on new XNP
      CALL D03PZW(XNP,CHISTR,NOP,XOP,M,FMON,NIP,I1,I1,IFAIL)
      IF (IFAIL.NE.0) GO TO 580
C     Max. integral of new monitor over each mesh interval
C     FMON on XNP held in CHISTR
      RMXINT = 0.0D0
      IDRMAX = 1
      DO 540 I = 2, NOP
         RINT = 0.5D0*(CHISTR(I-1)+CHISTR(I))*DXW(I-1)
         IF (RINT.LE.RMXINT) GO TO 540
         IDRMAX = I
         RMXINT = RINT
  540 CONTINUE
      WRITE (REC,FMT=99980)
      CALL X04BAF(IDEV,REC)
      WRITE (REC,FMT=99979) RMXINT, IDRMAX
      CALL X04BAF(IDEV,REC)
  560 IFDUM = IFAIL
      RETURN
C ----------------------------------------------------------------------
C     Error condition encountered
C ----------------------------------------------------------------------
  580 IFAIL = ISTAGE
      WRITE (REC,FMT=99984) ISTAGE
      CALL X04ABF(0,IDEV)
      CALL X04BAF(IDEV,REC)
C VP NOV94
C     If soft fail - returns to calling program
C     IF (IFDUM.NE.0) GO TO 560
C     Hard fail - prints useful info and then stops
      IF (ISTAGE.EQ.0 .OR. ISTAGE.EQ.4) IFDUM = -2
C END OF VP CHANGES
      IF (ISTAGE.EQ.0) THEN
         IF (IPMINF.GT.0) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAF(IDEV,REC)
            CALL D02NNN(XOP,NIP,29)
            WRITE (REC,FMT=99983)
            CALL X04BAF(IDEV,REC)
            CALL D02NNN(FMON,NIP,29)
         END IF
C VP NOV94 (ERROR MESSAGE CHANGED)
         ERRMSG =
     *' Invalid monitor function (FMON negative at one
     *  or more points) or invalid input mesh
     *  (XOP(NIP) not monotonically increasing). '
         CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
         RETURN
      ELSE
         WRITE (REC,FMT=99995)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(FMON,NIP,29)
      END IF
      IF (ISTAGE.LT.4) RETURN
C     ISTAGE=4 failure from inverse interpolation
      WRITE (REC,FMT=99986) NOP
      CALL X04BAF(IDEV,REC)
      CALL D02NNN(XNP,NOP,29)
      IF (ISTAGE.GT.4) THEN
         WRITE (REC,FMT=99982)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(CHISTR,NOP,29)
      END IF
      RETURN
C
99999 FORMAT (' **REMESH**  No. output points exceeds max. allowed ',I8,
     *       /'  sets nop=maxnpt and continues ')
99998 FORMAT (' INITIAL MESH')
99997 FORMAT (' WARNING: monitor function is zero everywhere - the mes',
     *       'h is unchanged')
99996 FORMAT (' REMESH CALLED AT TIME ',D12.4,' SCALING FACTOR FOR MON',
     *       'ITOR IS ',D12.4)
99995 FORMAT (' SCALED MONITOR FUNCTION')
99994 FORMAT (' PADDED MONITOR AFTER FORWARD SWEEP')
99993 FORMAT (' PADDED MONITOR AFTER REVERSE SWEEP')
99992 FORMAT (' MONITOR FUNCTION ZERO EVERYWHERE ')
99991 FORMAT (' NO. OF INPUT AND OUTPUT POINTS = ',2I4)
99990 FORMAT (' SUM AND CONST = ',2D12.4)
99989 FORMAT (' MAX ADJACENT MESH SIZE AT POINT',I5,' IS ',D12.4)
99988 FORMAT (' RATIO MAX:EVEN = ',D12.4)
99987 FORMAT (' OLD MESH')
99986 FORMAT (' NEW NO. PTS = ',I3,' NEW MESH =')
99985 FORMAT (' MESH SIZES ARE ')
99984 FORMAT (' *** ERROR IN REMESH ROUTINE AT STAGE ',I2,'  ***')
99983 FORMAT (' MONITOR FUNCTION')
99982 FORMAT (' NEW MONITOR ')
99981 FORMAT (' MAX. INTEGRAL OF SCALED MONITOR OVER EACH INPUT MESH I',
     *       'NTERVAL IS ')
99980 FORMAT (' MAX. INTEGRAL OF NEW MONITOR OVER EACH NEW MESH INTERV',
     *       'AL IS ')
99979 FORMAT (1X,D12.4,' AT POINT',I5)
      END
