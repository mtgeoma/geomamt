      SUBROUTINE D03PDR(PDEFCN,DUMPD1,NP,XP,UP,ITYPE,U,NPTS,NPDE,NEL,
     *                  NPTL,OMEGA,COEFF,XBK,IBK,IFLAG,NV,V,VDOT,RT,T,
     *                  IR,SPDEFN)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C ----------------------------------------------------------------------
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                        PARAMETER LIST
C                       ----------------
C
C     XP(NP)            The mesh points at which interpolated values
C                       are required. these points such be in
C                       increasing order.
C
C     UP(NPDE,NP,ITYPE)  Array that holds the values found by
C                        interpolation.
C
C     IF ITYPE >= 1      UP(J,K,1) holds the solution value at mesh
C                        point XP(K) for Jth PDE
C
C     IF ITYPE >= 2      UP(J,K,2) Holds the space deriv of the solution
C                        at point XP(K) for Jth PDE.
C
C     IF ITYPE >= 3      UP(J,K,3) holds the flux R at the point
C                        XP(K) for the Jth PDE.
C
C     U(NPDE,NPTS)       Original solution vector from the ODE code.
C
C     NPTS               The number of mesh points used in computing U.
C
C     NPDE               The number of PDEs in the problem.
C
C     NEL                The number of spatial elements in the mesh.
C
C     NPTL               The number of mesh points per element.
C                        Therefore NPTS = NEL*(NPTL-1) + 1.
C
C     OMEGA              Matrix used in mapping from the solution on a
C                        Spatial interval to its Chebyshev coeffs.
C
C     COEFFS             Workspace used to hold these coeffs.
C
C     XBK(IBK)           Array used to hold the breakpoints between the
C                        spatial elements.
C
C     IFLAG              Error flag set to 0 unless extrapolation is
C                        tried and then set to 1.
C
C     NV                 The size of the additional ode system that is
C                        coupled to the pde system.
C
C     V(NV)              Coupled ODE variables
C
C     VDOT(NV)           Time derivatives
C
C     T                  The current value of the time variable.
C
C NOTE:
C     These last four variables are only used if ITYPE = 3
C     otherwise dummy variables may be passed across.
C
C     IR                 IRES param to test for illegal values
C                        the method used is decomposition of
C                        the solution per element into chebyshev
C                        coefficients. This is done by matrix
C                        multiplication using the omega matrix .
C                        FFT could also be used.
C                        Interpolation is used to provide those
C                        solution values in the element (using
C                        Clenshaws algorithm).
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ---------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IBK, IFLAG, IR, ITYPE, NEL, NP, NPDE, NPTL,
     *                  NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  COEFF(NPDE,NPTL,2), OMEGA(NPTL,NPTL),
     *                  RT(NPDE,NPTL,3), U(NPDE,NPTS), UP(NPDE,NP,*),
     *                  V(*), VDOT(*), XBK(IBK), XP(NP)
C     .. Subroutine Arguments ..
      EXTERNAL          DUMPD1, PDEFCN, SPDEFN
C     .. Scalars in Common ..
      DOUBLE PRECISION  TWOU
      INTEGER           IIFLAG
C     .. Local Scalars ..
      DOUBLE PRECISION  AL, BR, BR1, BR2, TEM, TEM1
      INTEGER           I, II, IONE, IP, IP1, IX, IY, IZ, J, K, NM1
C     .. Local Arrays ..
      DOUBLE PRECISION  XCON(2)
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AD03PD/TWOU
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /AD03PD/
C     .. Executable Statements ..
C
C ... Treat each element separately ...
C
      TEM = 1.0D0 + TWOU
      TEM1 = 1.0D0 - TWOU
      IONE = 1
      IP = 0
      NM1 = NPTL - 1
      IZ = 0
      DO 280 I = 1, NEL
         IP1 = I + 1
   20    CONTINUE
         IP = IP + 1
         IF (IP.EQ.(NP+1)) GO TO 300
         IF (XP(IP).LT.(XBK(I)*TEM1-TWOU)) GO TO 20
         IF (XP(IP).GT.(XBK(I+1)*TEM+TWOU)) THEN
            IP = IP - 1
            GO TO 280
         END IF
C
         IF (XP(IP).GT.(XBK(I+1)*TEM1-TWOU)) THEN
            IF (I.LT.NEL .AND. ITYPE.GE.2) IZ = 1
         END IF
C
C     ------------------------------------------------
C ... Process a sequence of XP(J) values in element I ...
C ... IX = start of correct part of solution vector U ...
C ... form the Chebyshev coeffs in the array coeff    ...
C     ------------------------------------------------
C
         IX = NM1*(I-1)
         DO 80 K = 1, NPDE
            DO 60 J = 1, NPTL
               COEFF(K,J,1) = 0.0D0
               DO 40 II = 1, NPTL
                  COEFF(K,J,1) = COEFF(K,J,1) + OMEGA(J,II)*U(K,IX+II)
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
C
C ... Form the Chebyshev coeffs of the space deriv ...
C
         IF (ITYPE.GE.2) THEN
            DO 120 K = 1, NPDE
               COEFF(K,NPTL,2) = 0.0D0
               COEFF(K,NPTL-1,2) = 2.0D0*NM1*COEFF(K,NPTL,1)
               DO 100 J = 2, NM1
                  COEFF(K,NPTL-J,2) = COEFF(K,NPTL-J+2,2) + COEFF(K,
     *                                NPTL-J+1,1)*2*(NPTL-J)
  100          CONTINUE
               COEFF(K,1,2) = COEFF(K,1,2)*0.5D0
  120       CONTINUE
         END IF
C
         XCON(1) = 2.0D0/(XBK(I+1)-XBK(I))
         XCON(2) = -0.5D0*XCON(1)*(XBK(I+1)+XBK(I))
         IY = MIN(2,ITYPE)
  140    CONTINUE
         DO 200 II = 1, IY
            DO 180 K = 1, NPDE
               BR1 = 0.0D0
               BR2 = 0.0D0
C
C .... COEFF(K,NPTL) is the NPTLth  coeff of solution of PDE ...
C
               AL = (XP(IP)*XCON(1)+XCON(2))*2.0D0
               BR = COEFF(K,NPTL,II)
               DO 160 J = 1, NM1
                  BR2 = COEFF(K,NPTL-J,II) + AL*BR - BR1
                  BR1 = BR
                  BR = BR2
  160          CONTINUE
               IF (II.EQ.1) THEN
                  UP(K,IP,II) = BR - BR1*AL*0.5D0
C
               ELSE IF (IZ.LT.2) THEN
                  UP(K,IP,II) = (BR-BR1*AL*0.5D0)*XCON(1)
C
               ELSE
                  UP(K,IP,II) = 1.D0/(XBK(I+1)-XBK(I-1))*(UP(K,IP,II)
     *                          *(XBK(I)-XBK(I-1))+(BR-BR1*AL*0.5D0)
     *                          *XCON(1)*(XBK(I+1)-XBK(I)))
               END IF
C
  180       CONTINUE
  200    CONTINUE
C
C ... If required form the flux at the interploated points ...
C ... (unless deriv is being formed by weighted average in ...
C ... which case wait until the formation is complete.)    ...
C
         IF (ITYPE.GE.3 .AND. IZ.NE.1) THEN
C
C ... Zero workspaces used in the flux call ...
C
            DO 240 J = 1, 3
               DO 220 K = 1, NPDE
                  RT(K,1,J) = 0.0D0
  220          CONTINUE
  240       CONTINUE
            IR = 1
C
C ... Form the flux at the interpolated points ...
C
            CALL SPDEFN(PDEFCN,DUMPD1,T,XP(IP),IONE,NPDE,UP(1,IP,1),
     *                  UP(1,IP,2),RT(1,1,3),RT(1,1,2),UP(1,IP,3),NV,V,
     *                  VDOT,IR,IIFLAG)
C
            IF (IIFLAG.EQ.2) RETURN
C
C ... RT dimensions must be large enough here ...
C
            IF (IR.NE.1) GO TO 300
C
         END IF
C
         IF (IP.EQ.NP) GO TO 280
         IP = IP + 1
         IF (IZ.EQ.1) THEN
            IZ = 2
            GO TO 260
C
C ... To calculate the other elements contribution to deriv ...
C
         END IF
C
         IF (IZ.EQ.2) IZ = 0
         IF (XP(IP).LT.(XBK(I+1)*TEM1-TWOU)) THEN
C
C ... Process another point in this element ...
C
            GO TO 140
C
         ELSE IF (XP(IP).LT.(XBK(I+1)*TEM+TWOU)) THEN
            IF ((I+1).LT.NEL .AND. ITYPE.GE.2) IZ = 1
C
C ... IZ = 1 means that weighted average must be used for ...
C ... derivative values that are requested at XBK(I+1).   ...
C
            GO TO 140
C
         END IF
C
  260    CONTINUE
         IP = IP - 2
  280 CONTINUE
      RETURN
C
  300 CONTINUE
      IFLAG = 1
      RETURN
C
      END
