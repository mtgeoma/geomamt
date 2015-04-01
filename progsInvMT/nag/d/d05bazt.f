      SUBROUTINE D05BAZ(CG,CF,CK,ILIM1,ILIM2,ALIM,NSTART,NXTPOL,L1WT,
     *                  LSTWT,TOL,THRESH,H,WT,STWT,WKY1,WKY2,V,VF,VK,VG,
     *                  WVG,A,B,THETA,IQ,WWVG,WVK,NOUT,N,INFO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C       This subroutine starts computing and extrapolating  from
C       ILIM1 to ILIM2. If the required TOL is not satisfied then
C       ILIM1 := 2*ILIM1  and  ILIM2 := 2*ILIM2 and the  control
C       is given to the subroutine D05BAV.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALIM, H, THRESH, TOL
      INTEGER           ILIM1, ILIM2, INFO, IQ, L1WT, LSTWT, N, NOUT,
     *                  NSTART, NXTPOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(6,0:5), B(0:IQ-1), STWT(LSTWT,0:NSTART),
     *                  THETA(0:IQ-1), V(0:NSTART), VF(0:N), VG(0:N),
     *                  VK(0:N), WKY1(0:N/2), WKY2(0:N), WT(0:L1WT),
     *                  WVG(0:N/2), WVK(0:NSTART-1,0:IQ-1,0:IQ-1),
     *                  WWVG(0:NSTART,0:IQ-1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CF, CG, CK
      EXTERNAL          CF, CG, CK
C     .. Scalars in Common ..
      DOUBLE PRECISION  VGC
C     .. Local Scalars ..
      DOUBLE PRECISION  HXT, HXT2
      INTEGER           ILIM11, ILIM22, ILIMP1, INDEX, IUP, NDIM
C     .. External Subroutines ..
      EXTERNAL          D05BAL, D05BAQ, D05BAR, D05BAS, D05BAV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /BD05BA/VGC
C     .. Executable Statements ..
C
      HXT2 = H/NXTPOL
      HXT = 2.D0*HXT2
      ILIM11 = ILIM1/2
      ILIM22 = ILIM2/2
      VG(0) = HXT2*VGC
C
      CALL D05BAQ(WKY1,STWT,VK,WVG,VF,V,LSTWT,NSTART,ILIM11+1,ILIM22,2)
C
      CALL D05BAQ(WKY2,STWT,VK,VG,VF,V,LSTWT,NSTART,ILIM1+1,ILIM2,1)
C
      IF (ILIM11.GT.NSTART) CALL D05BAL(ILIM11+1,ILIM22,NSTART,WVG,VK,
     *                           WKY1,WT,L1WT,2)
C
      CALL D05BAL(ILIM1+1,ILIM2,NSTART,VG,VK,WKY2,WT,L1WT,1)
C
C     ---- Find the solution with  stepsize  H /2**(NXTPOL-1),  ----
C     ---- up to  H*NSTART.                                     ----
C
      CALL D05BAS(CG,HXT,ALIM,WT,L1WT,VK,WVG,WKY1,ILIM11,ILIM22,NSTART,
     *            2,INFO)
      IF (INFO.EQ.4) RETURN
C
C     ---- Find the solution with stepsize  H /2**NXTPOL, ----
C     ---- up to  H*NSTART.                               ----
C
      CALL D05BAS(CG,HXT2,ALIM,WT,L1WT,VK,VG,WKY2,ILIM1,ILIM2,NSTART,1,
     *            INFO)
C
C     ---- Check the errors  ----
C
      NDIM = ILIM22 - ILIM11
      ILIMP1 = ILIM11 + 1
C
      CALL D05BAR(WKY1(ILIMP1),WKY2(2*ILIM11+1),THRESH,NDIM,INDEX)
      IF (ABS(WKY1(ILIM11+INDEX)).GT.TOL) THEN
         IUP = ILIM2
         CALL D05BAV(CG,CF,CK,NXTPOL,VK,VG,WVG,VF,WKY1,WKY2,STWT,WT,V,A,
     *               B,THETA,IQ,WWVG,WVK,LSTWT,L1WT,NSTART,H,IUP,ALIM,
     *               TOL,THRESH,NOUT,N,INFO)
      END IF
      RETURN
      END
