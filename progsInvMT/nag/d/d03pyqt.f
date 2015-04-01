      SUBROUTINE D03PYQ(NP,XP,UP,ITYPE,U,NPTS,NPDE,NEL,NPTL,OMEGA,COEFF,
     *                  XBK,IBK,IFAIL1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C  ---------------------------------------------------------------------
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C         PARAMETER LIST
C         --------------
C         XP(NP)         The mesh points at which interpolated values
C                        are required. These points such be in
C                        increasing order.
C         UP(NPDE,NP,ITYPE)  Array that holds the values found by
C                            interpolation.
C         IF ITYPE >= 1  UP(j,k,1) Holds the solution value at mesh
C                                  point XP(k) FOR jth PDE.
C         IF ITYPE >= 2  UP(j,k,2) Holds the space deriv of the solution
C                                  at point XP(k) for jth PDE.
C
C         U(1..NEQN)     original solution vector from the ODE code.
C
C         NPTS           The number of mesh points used in computing U.
C         NPDE           The number of PDEs in the problem.
C         NEL            The number of spatial elements in the mesh.
C         NPTL           The number of mesh points per element.
C                        Therefore NPTS = NEL*(NPTL-1) + 1.
C         OMEGA          matrix used in mapping from the solution on a
C                        spatial interval to its Chebyshev coeffs.
C         COEFFS         Workspace used to hold these coeffs.
C         XBK(IBK)       Array used to hold the breakpoints between the
C                        spatial elements.
C         IFAIL1          Error flag set to 0 unless extrapolation is
C                        tried and then set to 1.
C
C         The method used is decomposition of the solution
C         per element into Chebyshev coefficients. This is done by
C         matrix multiplication using the omega matrix. FFT techniques
C         could also be used. Interpolation is used to provide those
C         Solution values in the element (using Clenshaws algorithm).
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ---------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IBK, IFAIL1, ITYPE, NEL, NP, NPDE, NPTL, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  COEFF(NPDE,NPTL,2), OMEGA(NPTL,NPTL),
     *                  U(NPDE,NPTS), UP(NPDE,NP,1), XBK(IBK), XP(NP)
C     .. Local Scalars ..
      DOUBLE PRECISION  AL, BR, BR1, BR2, CU, TEM, TEM1
      INTEGER           I, II, IP, IP1, IX, IY, IZ, J, K, NM1
C     .. Local Arrays ..
      DOUBLE PRECISION  XCON(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
      CU = X02AJF()
      TEM1 = 1 - CU
      TEM = 1 + CU
      IP = 0
      NM1 = NPTL - 1
      IZ = 0
      DO 240 I = 1, NEL
         IP1 = I + 1
   20    IP = IP + 1
C
         IF (IP.EQ.(NP+1)) GO TO 260
         IF (XP(IP).LT.(XBK(I)*TEM1-CU)) GO TO 20
         IF (XP(IP).GT.(XBK(I+1)*TEM+CU)) GO TO 220
         IF (XP(IP).GE.(XBK(I+1)*TEM1-CU)) THEN
            IF (I.LT.NEL .AND. ITYPE.GE.2) IZ = 1
         END IF
C
C        ... IZ = 1 means that weighted average must be used for ...
C        ... derivative values that are requested at XBK(I+1)    ...
C
C
C        ... Process a sequence of XP(j) values in element I ...
C        ... IX = start of correct part of solution vector U ...
C        ... form the Chebyshev coeffs in the array coeff.   ...
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
C        Form the Chebyshev coeffs of the space deriv.
C
         IF (ITYPE.EQ.2) THEN
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
  140    DO 200 II = 1, IY
            DO 180 K = 1, NPDE
               BR1 = 0.0D0
               BR2 = 0.0D0
C
C        ... COEFF(k,NPTL) is the NPTLth  coeff of solution of PDE k ...
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
         IF (IP.EQ.NP) GO TO 240
         IP = IP + 1
         IF (IZ.EQ.1) THEN
            IZ = 2
            GO TO 220
C
C           ... Calculate the other elements contribution to deriv. ...
C
         END IF
C
         IF (IZ.EQ.2) IZ = 0
C
         IF (XP(IP).LE.(XBK(I+1)*TEM1-CU)) THEN
            GO TO 140
         ELSE IF (XP(IP).LE.(XBK(I+1)*TEM+CU)) THEN
            IF (I.LT.NEL .AND. ITYPE.GE.2) IZ = 1
            GO TO 140
         END IF
C
C        ... IZ = 1 means that weighted average must be used for ...
C        ... derivative values that are requested at XBK(I+1)    ...
C
  220    IP = IP - 2
  240 CONTINUE
      RETURN
  260 IFAIL1 = 3
      RETURN
      END
