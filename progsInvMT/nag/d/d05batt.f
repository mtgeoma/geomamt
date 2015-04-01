      SUBROUTINE D05BAT(CG,ALIM,H,STWT,LSTWT,WKY,VF,VK,VG,INFO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     ---------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      This subroutine computes the starting value using the
C      Adam3 and BDF2 with a stepsize H. These methods only
C      require one starting value.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ---------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, H
      INTEGER           INFO, LSTWT
C     .. Array Arguments ..
      DOUBLE PRECISION  STWT(LSTWT,0:1), VF(0:1), VG(0:1), VK(0:1),
     *                  WKY(0:1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CG
      EXTERNAL          CG
C     .. Scalars in Common ..
      DOUBLE PRECISION  GTOL, VGC
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST, GEVAL, GFUN, HALIM, SCALE, VG0, WT11, Y0
      INTEGER           IFLAG, IND, IR
C     .. Local Arrays ..
      DOUBLE PRECISION  C(26)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          C05AXF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Common blocks ..
      COMMON            /AD05BA/GTOL
      COMMON            /BD05BA/VGC
C     .. Executable Statements ..
C
C
C     ---- Evaluation of starting values with H  ----
C
C
      WKY(0) = VF(0)
      VG0 = H*VGC
      VG(0) = VG0
      CONST = VF(1) + STWT(1,0)*VK(1)*VG0
      WT11 = H*STWT(1,1)*VK(0)
      SCALE = SQRT(X02AJF())
      IR = 1
      IFLAG = 1
      Y0 = WKY(0)
      IND = 1
   20 CALL C05AXF(Y0,GFUN,GTOL,IR,SCALE,C,IND,IFLAG)
      GEVAL = CG(HALIM,Y0)
      IF (IND.EQ.0) GO TO 40
      GFUN = Y0 - WT11*GEVAL - CONST
      GO TO 20
   40 CONTINUE
      IF (IFLAG.GT.0) THEN
         INFO = 4
         RETURN
      END IF
      WKY(1) = Y0
      VG(1) = H*GEVAL
      RETURN
      END
