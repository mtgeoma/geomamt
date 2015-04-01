      SUBROUTINE D03PHL(PHFPDE,PHFBND,PDEPCF,BNDPCF,X,T,U,UDOT,NPDE,
     *                  NPTS,NC,RES,BETA,GAMMA,UMEAN,UX,F,G,CL,CF,RFUN,
     *                  IRES,V,VDOT,NV,FT,BT,GT,VT,IR,SPDEF1,BNDR,APLUS,
     *                  BPLUS,EPLUS,XHAT,PR,DPL,DPR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   This subroutine serves as an interface between a
C   parabolic partial differential equation in one
C   space variable and an integrator of initial value problems,
C
C   The differential equation should be of the form
C
C     C(I,X,U,UX,V)*(D/DT)U(I) = (D/DX)(X**NC*F(I,X,T,U,UX,V)))/X**NC
C                                                  .
C                                 - G(I,X,T,U,UX,V,V),I= 1,...,NPDE;
C
C   with boundary conditions
C                                              .
C     BETA(I)*F(I,X,U,UX,V) = GAMMA(I,U,UX,V,V,V) ,I = 1,...,NPDE;
C                           AT X = X(1) OR X(NPTS)
C
C   D03PHL transforms the system of partial differential equations
C   into a system of ordinary equations by semi-discretization in
C   the space variable x; this semi-discretization is performed
C   by application of the Skeel and Berzins (1987) discretisation
C
C  INPUTS:
C
C   X(NPTS); X(1),...,X(NPTS) IS A PARTITION OF (X(1),X(NPTS))
C
C   T  ;  THE TIME VARIABLE;
C
C   U(NPDE,NPTS); U(I,J) is an approximation of U(I,X(J),T),
C                 I = 1,...,NPDE; J = 1,...,NPTS;
C
C   NPDE ;   the number of partial differential equations
C
C   NPTS ;   the number of gridpoints;
C
C   NC   ;  A number designing the kind of space coordinates;
C           NC = 0 : Cartesian coordinates;
C           NC = 1 : circular coordinates;
C           NC = 2 : spherical coordinates;
C
C   RES(NPDE,NPTS) ; empty array
C
C     BETA,GAMMA,UMEAN,UX,F,G,CT,FT,GT,BT
C     work arrays of dimension NPDE
C     VT : work array of dimension NV.
C
C   OUTPUT:
C
C   RES(NPDE,NPTS); the  residual of the semi-discretised PDE.
C
C   N.B. to handle IRES = -1 the following options are used.
C
C   IF IRES = -1 and NV = 0 then IT=0 and NVIT = 1 to calculate only
C   those parts of resid depending on the time deriv
C   in the case when NV non-zero then IT = 1 and NVIT only is zero
C   when IRES is -1 .
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RFUN, T
      INTEGER           IR, IRES, NC, NPDE, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  APLUS(NPTS), BETA(NPDE), BPLUS(NPTS),
     *                  BT(NPDE,NPDE), CF(NPDE,NPDE), CL(NPDE,NPDE),
     *                  DPL(NPTS), DPR(NPTS), EPLUS(NPTS), F(NPDE),
     *                  FT(NPDE), G(NPDE), GAMMA(NPDE), GT(NPDE),
     *                  PR(NPTS), RES(NPDE,NPTS), U(NPDE,NPTS),
     *                  UDOT(NPDE,NPTS), UMEAN(NPDE), UX(NPDE), V(*),
     *                  VDOT(*), VT(*), X(NPTS), XHAT(NPTS)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPCF, BNDR, PDEPCF, PHFBND, PHFPDE, SPDEF1
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           IIFLAG, IOVFLO
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, XPOW
      INTEGER           I, IBND, IT, J, JTT, K, L, NVIT
C     .. Local Arrays ..
      INTEGER           IZ(3)
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /CD03PC/
C     .. Executable Statements ..
C
      JTT = 1
      IF (NC.GT.0 .AND. X(1).LE.UROUND) JTT = 0
C
C ... JTT = 0  indicates polar singularity at X = 0 being treated ...
C
      IT = 1
      NVIT = 1
      IF (IRES.EQ.-1) THEN
         IT = 0
         IF (NV.GT.0) THEN
            NVIT = 0
            IT = 1
         END IF
      END IF
      TEMP = 0.0D0
      DO 20 I = 1, 3
         IZ(I) = 1
   20 CONTINUE
C
C ... First segment,process the boundary conditions; if of Neumann ...
C ... type the value of the flux F at the boundary is updated.     ...
C
      DO 40 I = 1, NPDE
         RES(I,1) = 0.0D0
         UX(I) = U(I,2)*DPR(2) + U(I,1)*DPL(2)
   40 CONTINUE
      IBND = 0
      CALL BNDR(PHFBND,BNDPCF,T,BETA,GAMMA,U,UX,NPDE,IBND,NV,V,VDOT,
     *          IZ(1),IIFLAG)
C
      IF (IIFLAG.EQ.2) THEN
         IRES = IZ(1)
         RETURN
      END IF
C
      IF (NVIT.EQ.0) THEN
C
         CALL BNDR(PHFBND,BNDPCF,T,BETA,BT,U,UX,NPDE,IBND,NV,V,VT,IZ(1),
     *             IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(1)
            RETURN
         END IF
C
         DO 60 I = 1, NPDE
            GAMMA(I) = GAMMA(I) - BT(I,1)
   60    CONTINUE
      END IF
      XPOW = 1.D0
      IF (NC.GT.0) XPOW = X(1)**NC
      DO 100 J = 1, NPDE
         IF (BETA(J).NE.0.0D0) RES(J,1) = -XPOW*GAMMA(J)/BETA(J)*IT
         DO 80 K = 1, NPDE
            CL(J,K) = 0.0D0
   80    CONTINUE
  100 CONTINUE
C
C ... Here the implementation of the left boundary conditions ...
C ... ends and  the real assembly begins.                     ...
C
      DO 240 L = 2, NPTS
         DO 140 J = 1, NPDE
            UMEAN(J) = PR(L)*(U(J,L)-U(J,L-1)) + U(J,L-1)
            UX(J) = DPL(L)*U(J,L-1) + DPR(L)*U(J,L)
            DO 120 K = 1, NPDE
               CF(K,J) = 0.0D0
  120       CONTINUE
  140    CONTINUE
         CALL SPDEF1(PHFPDE,PDEPCF,T,XHAT(L),NPDE,UMEAN,UX,CF,G,F,NV,V,
     *               VDOT,IZ(2),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(2)
            RETURN
         END IF
C
         IF (NVIT.EQ.0) THEN
            CALL SPDEF1(PHFPDE,PDEPCF,T,XHAT(L),NPDE,UMEAN,UX,BT,GT,FT,
     *                  NV,V,VT,IZ(2),IIFLAG)
C
            IF (IIFLAG.EQ.2) THEN
               IRES = IZ(2)
               RETURN
            END IF
C
            DO 180 I = 1, NPDE
               G(I) = G(I) - GT(I)
               F(I) = F(I) - FT(I)
               TEMP = TEMP + F(I)
               DO 160 J = 1, NPDE
                  TEMP = TEMP + (CF(J,I)-BT(J,I))
  160          CONTINUE
  180       CONTINUE
         END IF
         DO 220 J = 1, NPDE
            RES(J,L-1) = (RES(J,L-1)-G(J)*APLUS(L)+F(J)*EPLUS(L))
     *                   /(APLUS(L)+BPLUS(L-1))*IT
            RES(J,L) = -G(J)*BPLUS(L) - F(J)*EPLUS(L)*JTT
            DO 200 K = 1, NPDE
               RES(J,L-1) = RES(J,L-1) - (APLUS(L)*CF(J,K)+CL(J,K))
     *                      /(APLUS(L)+BPLUS(L-1))*UDOT(K,L-1)
               CL(J,K) = BPLUS(L)*CF(J,K)
  200       CONTINUE
  220    CONTINUE
         JTT = 1
  240 CONTINUE
C
C ... Finish processing any non-flux boundary conditions at X(1) ...
C
      DO 260 J = 1, NPDE
         IF (BETA(J).EQ.0.D0) RES(J,1) = GAMMA(J)*IT
  260 CONTINUE
C
C ... Finally, process the right boundary conditions ...
C
      XPOW = 1.0D0
      IF (NC.GT.0) XPOW = X(NPTS)**NC
      IBND = 1
      CALL BNDR(PHFBND,BNDPCF,T,BETA,GAMMA,U(1,NPTS),UX,NPDE,IBND,NV,V,
     *          VDOT,IZ(3),IIFLAG)
C
      IF (IIFLAG.EQ.2) THEN
         IRES = IZ(3)
         RETURN
      END IF
C
      IF (NVIT.EQ.0) THEN
         CALL BNDR(PHFBND,BNDPCF,T,BETA,BT,U(1,NPTS),UX,NPDE,IBND,NV,V,
     *             VT,IZ(3),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(3)
            RETURN
         END IF
C
         DO 280 I = 1, NPDE
            GAMMA(I) = GAMMA(I) - BT(I,1)
  280    CONTINUE
      END IF
      DO 320 J = 1, NPDE
         IF (BETA(J).EQ.0.D0) THEN
            RES(J,NPTS) = GAMMA(J)*IT
         ELSE
            RES(J,NPTS) = (RES(J,NPTS)+XPOW*GAMMA(J)/BETA(J))
     *                    *IT/BPLUS(NPTS)
            DO 300 K = 1, NPDE
               RES(J,NPTS) = RES(J,NPTS) - CF(J,K)*UDOT(K,NPTS)
  300       CONTINUE
         END IF
  320 CONTINUE
C
C ... Illegal problem spec encountered ...
C
      IF (ABS(TEMP).GE.UROUND) THEN
         IIFLAG = 1
         RETURN
      END IF
      DO 340 I = 1, 3
         IF (IZ(I).NE.1) THEN
            IR = IZ(I)
            GO TO 360
         END IF
  340 CONTINUE
  360 CONTINUE
      RETURN
      END
