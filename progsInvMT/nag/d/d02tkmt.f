      SUBROUTINE D02TKM(XCOL,HRHO,JCOL,KCOL,WI,VI,IPVTW,KD,ZVAL,DF,ACOL,
     *                  DMZO,FJAC,MSING,M,NEQ,MAXORD,MSTAR,SETRES)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   Purpose:
C      Construct a group of  NEQ  rows of the matrices  WI  and  VI.
C      corresponding to an interior collocation point of the i-th
C      subinterval
C
C   Arguments:
C      XCOL   - the location of the collocation point
C      HRHO   - length of i-th subinterval times relative fraction
C               of collocation point along interval (ie. XCOL - x_i)
C      JCOL   - XCOL is the JCOL-th collocation point
C               in the i-th subinterval.
C      KCOL   - the number of collocation points in the i-th
C               subinterval.
C      WI,VI  - part of the i-th block of the collocation matrix
C               before parameter condensation.
C               WI is (NEQ*KCOL) x (NEQ*KCOL),
C               VI is (NEQ*KCOL) x (MSTAR)
C      IPVTW  - pivot info for decomposition of WI (required later)
C      KD     - no. of rows in VI and WI . (NEQ*KCOL)
C      ZVAL   - z(XCOL) of length MSTAR but stored here as a 2-d array
C               of dimension (NEQ,*)
C      DF     - the jacobian at XCOL .
C      ACOL   - the mesh independent rk-basis coefficients
C      DMZO   - the residual of the ode system
c      FJAC  - external to evaluate the Jacobian
C      MSING  - flag for singularity in WI
C      M      - the order of each equation
C      NEQ    - the numder of odes
C      MAXORD - the maximal order of the odes
C      SETRES - determine if residual od system should be computed
C               (for a new mesh at the first iteration of a nonlinear
C               system)
C
C   Author:
C      R.W. Brankin, NAG Ltd., Aug 1994
C      (modified from COLNEW routine VWBLOK)
C
C**********************************************************************
C
C For the first collocation point, initialize WI to indentity
C Note that WI is the 0 matrix on entry in this instance
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HRHO, XCOL
      INTEGER           JCOL, KCOL, KD, MAXORD, MSING, MSTAR, NEQ
      LOGICAL           SETRES
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOL(7,4), DF(NEQ,NEQ,MAXORD), DMZO(KD),
     *                  VI(KD,MSTAR), WI(KD,KD), ZVAL(NEQ,MAXORD)
      INTEGER           IPVTW(KD), M(NEQ)
C     .. Subroutine Arguments ..
      EXTERNAL          FJAC
C     .. Local Scalars ..
      DOUBLE PRECISION  AJL, BL, FACT
      INTEGER           I, I0, I1, I2, ID, IR, IW, J, J1, JCOMP, JDF,
     *                  JN, JV, JW, L, LL, LP1, MJ
C     .. Local Arrays ..
      DOUBLE PRECISION  BASM(5), HA(7,4)
C     .. External Subroutines ..
      EXTERNAL          DGETRF, DGETRS
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IF (JCOL.EQ.1) THEN
         DO 20 ID = 1, KD
            WI(ID,ID) = 1.D0
   20    CONTINUE
      END IF
C
C Calculate local basis
C
      FACT = 1.D0
      DO 60 L = 1, MAXORD
         FACT = FACT*HRHO/DBLE(L)
         BASM(L) = FACT
         DO 40 J = 1, KCOL
            HA(J,L) = FACT*ACOL(J,L)
   40    CONTINUE
   60 CONTINUE
C
C Zero jacobian prior to evaluation
C
      DO 120 J1 = 1, MAXORD
         DO 100 IR = 1, NEQ
            DO 80 J = 1, NEQ
               DF(J,IR,J1) = 0.D0
   80       CONTINUE
  100    CONTINUE
  120 CONTINUE
C
      CALL FJAC(XCOL,ZVAL,NEQ,M,DF)
C
C  Build neq rows for each interior collocation point x.
C  The linear expressions to be constructed are:
C   (m(id))
C  u     -  df(id,1)*z(1) - ... - df(id,mstar)*z(mstar)
C   id
C  for id = 1 to neq.
C
      I0 = (JCOL-1)*NEQ
      I1 = I0 + 1
      I2 = I0 + NEQ
C
C Evaluate  dmzo = dmz - df * zval  once for a new mesh
C
      IF (SETRES) THEN
         DO 180 I = 1, NEQ
            DO 160 J = 1, M(I)
               FACT = -ZVAL(I,J)
               DO 140 ID = 1, NEQ
                  DMZO(I0+ID) = DMZO(I0+ID) + FACT*DF(ID,I,J)
  140          CONTINUE
  160       CONTINUE
  180    CONTINUE
      END IF
C
C Loop over the  neq  expressions to be set up for the
C current collocation point.
C
      J1 = 0
      DO 240 I = 1, NEQ
         DO 220 J = 1, M(I)
            J1 = J1 + 1
            DO 200 ID = 1, NEQ
               VI(I0+ID,J1) = DF(ID,I,J)
  200       CONTINUE
  220    CONTINUE
  240 CONTINUE
      JN = 1
      DO 360 JCOMP = 1, NEQ
         MJ = M(JCOMP)
         JN = JN + MJ
         DO 340 L = 1, MJ
            JV = JN - L
            JW = JCOMP
            DO 280 J = 1, KCOL
               AJL = -HA(J,L)
               DO 260 IW = I1, I2
                  WI(IW,JW) = WI(IW,JW) + AJL*VI(IW,JV)
  260          CONTINUE
               JW = JW + NEQ
  280       CONTINUE
            LP1 = L + 1
            IF (L.LT.MJ) THEN
               DO 320 LL = LP1, MJ
                  JDF = JN - LL
                  BL = BASM(LL-L)
                  DO 300 IW = I1, I2
                     VI(IW,JV) = VI(IW,JV) + BL*VI(IW,JDF)
  300             CONTINUE
  320          CONTINUE
            END IF
  340    CONTINUE
  360 CONTINUE
C
C If all the collocation points for this interval have been processed
C then do parameter condensation, that is decompose the WI block
C and solve for the MSTAR columns of VI
C
      IF (JCOL.EQ.KCOL) THEN
         MSING = 0
         CALL DGETRF(KD,KD,WI,KD,IPVTW,MSING)
         IF (MSING.EQ.0) THEN
            CALL DGETRS('N',KD,MSTAR,WI,KD,IPVTW,VI,KD,MSING)
         END IF
      END IF
      RETURN
      END
