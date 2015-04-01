      SUBROUTINE D03PDS(PDEFCN,BNDARY,DUMPD1,DUMPD2,NPDE,NPTS,T,U,RES,
     *                  UDOT,M,X,OMEGA,DU,XBK,BETA,GAMMA,DUDX,R,Q,NEL,
     *                  NPTL,XC,CCR,IRES,RT,QT,UDT,C,NV,V,VDOT,VDUM,
     *                  SPDEFN,SBNDR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C    -------------------------------------------------------------------
C    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C       CHEBYSHEV C0 COLLOCATION ROUTINE
C       REFERENCE  ; C0 CHEBYSHEV METHODS FOR PARABOLIC EQUATIONS.
C                    BY M. BERZINS AND P.M. DEW
C    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    -------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, M, NEL, NPDE, NPTL, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  BETA(NPDE,4), C(NPDE,NPDE,NPTL), CCR(NPTL),
     *                  DU(NPTL,NPTL), DUDX(NPDE,NPTL), GAMMA(NPDE,4),
     *                  OMEGA(NPTL,NPTL), Q(NPDE,NPTL), QT(NPDE,NPTL),
     *                  R(NPDE,NPTL), RES(NPDE,NPTS), RT(NPDE,NPTL),
     *                  U(NPDE,NPTS), UDOT(NPDE,NPTS), UDT(NPDE,NPTL),
     *                  V(*), VDOT(*), VDUM(*), X(NPTS), XBK(*),
     *                  XC(NPTL)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, DUMPD1, DUMPD2, PDEFCN, SBNDR, SPDEFN
C     .. Scalars in Common ..
      DOUBLE PRECISION  TWOU
      INTEGER           IDEV, IIFLAG, ITRACE
C     .. Local Scalars ..
      DOUBLE PRECISION  CCLF, CCRT, H, MP1, SAVEL, SAVER, SFIRST, SUM,
     *                  TEM
      INTEGER           I, II, IJ, IK, IT, IV, J, JJ, JK, K, KJ, NM1
C     .. Local Arrays ..
      INTEGER           IZ(6)
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /AD03PD/TWOU
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /AD02NM/, /AD03PD/
C     .. Executable Statements ..
      NM1 = NPTL - 1
      IT = 1
      IV = MAX(1,NV)
      SUM = 0.0D0
      MP1 = 1.0D0
      DO 460 I = 1, NEL
         JJ = (I-1)*NM1
         IJ = JJ + 1
         H = 2.0D0/(XBK(I+1)-XBK(I))
         DO 20 IK = 1, 6
            IZ(IK) = 1
   20    CONTINUE
C
C ... Main loop over all the spatial elements start by forming the ...
C ... space derivs of U and UDOT in DUDX and UTDX respectively.    ...
C
         DO 80 K = 1, NPDE
            DO 60 II = 1, NPTL
               DUDX(K,II) = 0.0D0
               DO 40 J = 1, NPTL
                  DUDX(K,II) = DUDX(K,II) + DU(II,J)*U(K,JJ+J)*H
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
C
         IF (I.EQ.1) THEN
C
C ... Save the values needed for left boundary conditions ...
C
            DO 100 J = 1, NPDE
               BETA(J,3) = DUDX(J,1)
  100       CONTINUE
         END IF
C
         IF (I.EQ.NEL) THEN
C
C ... Save the values needed for right boundary conditions ...
C
            DO 120 J = 1, NPDE
               GAMMA(J,3) = DUDX(J,NPTL)
  120       CONTINUE
         END IF
C
C ... Evaluate the functions Q and R in this element ...
C
         CALL SPDEFN(PDEFCN,DUMPD1,T,X(IJ),NPTL,NPDE,U(1,IJ),DUDX,C,Q,R,
     *               IV,V,VDOT,IZ(1),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(1)
            RETURN
         END IF
C
         IF (IRES.EQ.-1) THEN
C
C ... Form those parts of Q and R involving UDOT ...
C
            DO 160 J = 1, NPDE
               DO 140 K = 1, NPTL
                  UDT(J,K) = 0.0D0
  140          CONTINUE
  160       CONTINUE
            CALL SPDEFN(PDEFCN,DUMPD1,T,X(IJ),NPTL,NPDE,U(1,IJ),DUDX,C,
     *                  QT,RT,IV,V,VDUM,IZ(2),IIFLAG)
C
            IF (IIFLAG.EQ.2) THEN
               IRES = IZ(2)
               RETURN
            END IF
C
            DO 200 J = 1, NPDE
               DO 180 K = 1, NPTL
                  Q(J,K) = Q(J,K) - QT(J,K)
                  R(J,K) = R(J,K) - RT(J,K)
                  SUM = SUM + ABS(R(J,K))
  180          CONTINUE
  200       CONTINUE
C
C ... Note if SUM not equal to zero, it means that flux ...
C ... function illegally depends on time derivatives    ...
C
            IF (SUM.GT.TWOU) THEN
               IIFLAG = 1
               GO TO 540
            END IF
         END IF
C
C ... Form proper Q function as in paper by adding C UDOT term ...
C
         DO 260 J = 1, NPTL
            DO 240 K = 1, NPDE
               DO 220 KJ = 1, NPDE
                  Q(K,J) = Q(K,J) + C(K,KJ,NPTL)*UDOT(KJ,IJ+J-1)
  220          CONTINUE
  240       CONTINUE
  260    CONTINUE
         IF (M.GT.0) THEN
C
C ... Modify Q function if polar co-ordinates ...
C
            KJ = 1
            IF (X(IJ).LE.TWOU) THEN
               MP1 = 1.0D0 + M
               KJ = 2
               DO 280 K = 1, NPDE
C                   R(K,1) = 0.0D0
                  Q(K,1) = Q(K,1)/(M+1)
  280          CONTINUE
            END IF
C
            DO 320 J = KJ, NPTL
               DO 300 K = 1, NPDE
                  Q(K,J) = Q(K,J) - R(K,J)*M/X(JJ+J)
  300          CONTINUE
  320       CONTINUE
         END IF
C
C ... Set up savel and saver  for the boundary and interface ...
C ... conditions and form drdx  by overwriting DUDX.         ...
C
         KJ = MAX(2,I)
         JK = MIN(NEL,I+1) + 1
         SAVEL = 1.0D0/(XBK(KJ)+XBK(I+1)-XBK(KJ-1)-XBK(I))
         SAVER = 1.0D0/(XBK(JK)+XBK(I+1)-XBK(JK-1)-XBK(I))
         IF (I.EQ.1) SFIRST = SAVEL
         DO 380 K = 1, NPDE
            DO 360 II = 2, NM1
               DUDX(K,II) = 0.0D0
               DO 340 J = 1, NPTL
                  DUDX(K,II) = DUDX(K,II) + DU(II,J)*R(K,J)
  340          CONTINUE
  360       CONTINUE
  380    CONTINUE
C
C ... Form the residual and the interface conditions ..
C
         DO 420 J = 1, NPDE
            CCLF = 0.0D0
            CCRT = 0.0D0
            DO 400 K = 2, NM1
               RES(J,JJ+K) = Q(J,K) - DUDX(J,K)*H
C
C ... Clenshaw-Curtis components for boundaries  ...
C
               CCLF = CCLF + CCR(K)*(DUDX(J,K)*(1.0D0-XC(K))-R(J,K))
               CCRT = CCRT + CCR(K)*(DUDX(J,K)*(XC(K)+1.0D0)+R(J,K))
  400       CONTINUE
            JK = IJ + NM1
            TEM = R(J,1) + R(J,NPTL)
            RES(J,IJ) = RES(J,IJ) + (Q(J,1)*2.0D0/H-TEM+CCLF)*SAVEL
            RES(J,JK) = (Q(J,NPTL)*2.0D0/H+TEM+CCRT)*SAVER
  420    CONTINUE
C
         DO 440 IK = 1, 2
            IF (IZ(IK).NE.1) THEN
               IRES = IZ(IK)
               GO TO 540
            END IF
C
  440    CONTINUE
  460 CONTINUE
C
C ... Evaluate the functions BETA and GAMMA at the ...
C ... boundary conditions                          ...
C
      CALL SBNDR(BNDARY,DUMPD2,T,BETA(1,1),GAMMA(1,1),U(1,1),BETA(1,3),
     *           NPDE,.TRUE.,IV,V,VDOT,IZ(3),IIFLAG)
C
      IF (IIFLAG.EQ.2) THEN
         IRES = IZ(3)
         RETURN
      END IF
C
      CALL SBNDR(BNDARY,DUMPD2,T,BETA(1,2),GAMMA(1,2),U(1,NPTS),
     *           GAMMA(1,3),NPDE,.FALSE.,IV,V,VDOT,IZ(4),IIFLAG)
C
      IF (IIFLAG.EQ.2) THEN
         IRES = IZ(4)
         RETURN
      END IF
C
      IF (IRES.EQ.-1) THEN
C
C ... Isolate parts of GAMMA(J,I) I = 1,2  depending on VDOT,YDOT ...
C
         CALL SBNDR(BNDARY,DUMPD2,T,BETA(1,1),BETA(1,4),U(1,1),BETA(1,3)
     *              ,NPDE,.TRUE.,IV,V,VDUM,IZ(5),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(5)
            RETURN
         END IF
C
         CALL SBNDR(BNDARY,DUMPD2,T,BETA(1,2),GAMMA(1,4),U(1,NPTS),
     *              GAMMA(1,3),NPDE,.FALSE.,IV,V,VDUM,IZ(6),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(6)
            RETURN
         END IF
C
         DO 480 J = 1, NPDE
            GAMMA(J,1) = GAMMA(J,1) - BETA(J,4)
            GAMMA(J,2) = GAMMA(J,2) - GAMMA(J,4)
  480    CONTINUE
      END IF
C
C ... LH and RH boundary conditions are processed ...
C
      DO 500 J = 1, NPDE
         RES(J,1) = MP1*(RES(J,1)*BETA(J,1)*2.0D0+GAMMA(J,1)
     *              *4.0D0/CCR(1)*SFIRST)
         RES(J,NPTS) = RES(J,NPTS)*BETA(J,2)*2.0D0 - GAMMA(J,2)
     *                 *4.0D0/CCR(1)*SAVER
  500 CONTINUE
      DO 520 IK = 3, 6
         IF (IZ(IK).NE.1) IRES = IZ(IK)
  520 CONTINUE
  540 CONTINUE
      RETURN
      END
