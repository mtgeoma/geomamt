      SUBROUTINE D03PLL(PLFPDE,PLFBND,PDEPFF,BNDPFF,X,T,U,UDOT,NPDE,
     *                  NPTS,RES,BETA,G,UMEAN,UX,F,S,PL,PF,IRES,V,VDOT,
     *                  NV,FT,BT,ST,DC,DD,UB,CL,CF,CS,DS,VT,IR,RFLUX,
     *                  SPDEF1,BNDR,APLUS,BPLUS,EPLUS,XHAT,PR,DPL,DPR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   This subroutine serves as an interface between a system of PDEs
C   of convection/diffusion type in one space variable and an integrator
C   of initial value problems.
C
C   The PDEs should be of the general form
C
C     P(i,j)*dU(j)/dt + dF(i)/dx = C(i)*dD(i)/dx + S(i)
C
C     where F = convective flux term,
C           C*dD/dx = diffusive term,
C           S = source term.
C
C   with boundary conditions
C               .
C     G(I,X,U,V,V), I = 1,...,NPDE; at X = X(1) or X(NPTS).
C
C   D03PLL transforms the system of partial differential equations
C   into a system of ordinary equations by semi-discretization in
C   the space variable x; this semi-discretization is performed
C   by application of the Skeel and Berzins (1987) discretisation
C   with upwind differencing of the convective flux terms.
C   Ref: Pennington and Berzins, ACM TOMS, vol 20, pp.63-99.
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
C   RES(NPDE,NPTS) ; empty array
C
C   RFLUX(NPDE,NPTS) ; array of convective flux values at midpoints
C
C     BETA,G,UMEAN,UX,F,S,FT,ST,DC,DD,CL,CF,CT,DT,CS,DS: work arrays
C     of dimension NPDE.
C     PF,PL,BT: work arrays of dimension NPDE*NPDE.
C     UB : work array of dimension 3*NPDE.
C     VT : work array of dimension NV.
C
C   OUTPUT:
C
C   RES(NPDE,NPTS); the  residual of the semi-discretised PDE.
C
C   N.B. to handle IRES = -1 the following options are used:
C   If IRES = -1 and NV = 0 then IT=0 and NVIT = 1 to calculate only
C   those parts of RESID depending on the time derivative.
C   In the case when NV non-zero then IT = 1 and NVIT only is zero
C   when IRES is -1 .
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IR, IRES, NPDE, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  APLUS(NPTS), BETA(NPDE), BPLUS(NPTS),
     *                  BT(NPDE,NPDE), CF(NPDE), CL(NPDE), CS(NPDE),
     *                  DC(NPDE), DD(NPDE), DPL(NPTS), DPR(NPTS),
     *                  DS(NPDE), EPLUS(NPTS), F(NPDE), FT(NPDE),
     *                  G(NPDE), PF(NPDE,NPDE), PL(NPDE,NPDE), PR(NPTS),
     *                  RES(NPDE,NPTS), RFLUX(NPDE,NPTS), S(NPDE),
     *                  ST(NPDE), U(NPDE,NPTS), UB(NPDE,3),
     *                  UDOT(NPDE,NPTS), UMEAN(NPDE), UX(NPDE), V(*),
     *                  VDOT(*), VT(*), X(NPTS), XHAT(NPTS)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPFF, BNDR, PDEPFF, PLFBND, PLFPDE, SPDEF1
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           IIFLAG, IOVFLO
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, IBND, IT, J, K, L, NVIT
C     .. Local Arrays ..
      INTEGER           IZ(3)
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /CD03PC/, /XD03PC/
C     .. Executable Statements ..
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
C ... First segment -- process the left-hand boundary conditions ...
C
      DO 40 J = 1, NPDE
         RES(J,1) = 0.0D0
   40 CONTINUE
      IBND = 0
      CALL BNDR(PLFBND,BNDPFF,T,G,U,UB,X,NPDE,NPTS,IBND,NV,V,VDOT,IZ(1),
     *          IIFLAG)
C
      IF (IIFLAG.EQ.2) THEN
         IRES = IZ(1)
         RETURN
      END IF
C
      IF (NVIT.EQ.0) THEN
C
         CALL BNDR(PLFBND,BNDPFF,T,BT,U,UB,X,NPDE,NPTS,IBND,NV,V,VT,
     *             IZ(1),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(1)
            RETURN
         END IF
C
         DO 60 J = 1, NPDE
            G(J) = G(J) - BT(J,1)
   60    CONTINUE
      END IF
C
      DO 100 J = 1, NPDE
         CS(J) = 0.0D0
         DS(J) = 0.0D0
         DO 80 K = 1, NPDE
            PL(J,K) = 0.0D0
   80    CONTINUE
  100 CONTINUE
C
C ... Here the implementation of the left boundary conditions ends and
C     the real assembly begins ...
C
      DO 240 L = 2, NPTS
         DO 140 J = 1, NPDE
            UMEAN(J) = PR(L)*(U(J,L)-U(J,L-1)) + U(J,L-1)
            UX(J) = DPL(L)*U(J,L-1) + DPR(L)*U(J,L)
            DO 120 K = 1, NPDE
               PF(K,J) = 0.0D0
  120       CONTINUE
  140    CONTINUE
         CALL SPDEF1(PLFPDE,PDEPFF,T,XHAT(L),NPDE,UMEAN,UX,PF,CF,DC,S,
     *               NV,V,VDOT,IZ(2),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(2)
            RETURN
         END IF
C
         IF (NVIT.EQ.0) THEN
            CALL SPDEF1(PLFPDE,PDEPFF,T,XHAT(L),NPDE,UMEAN,UX,BT,CL,DD,
     *                  ST,NV,V,VT,IZ(2),IIFLAG)
C
            IF (IIFLAG.EQ.2) THEN
               IRES = IZ(2)
               RETURN
            END IF
C
            DO 180 J = 1, NPDE
               S(J) = S(J) - ST(J)
               DC(J) = DC(J) - DD(J)
               CF(J) = CF(J) - CL(J)
               TEMP = TEMP + DC(J) + CF(J)
               DO 160 K = 1, NPDE
                  TEMP = TEMP + (PF(K,J)-BT(K,J))
  160          CONTINUE
  180       CONTINUE
         END IF
C
C ... Calculate residuals ...
C
         DO 220 J = 1, NPDE
C  NB Next expression valid only for cartesian coords (EPLUS=1).
            RES(J,L-1) = (RES(J,L-1)+S(J)*APLUS(L)-RFLUX(J,L)*EPLUS(L)
     *                   +(BPLUS(L-1)*CS(J)+APLUS(L)*CF(J))*(DC(J)-DS(J)
     *                   )/(APLUS(L)+BPLUS(L-1)))/(APLUS(L)+BPLUS(L-1))
     *                   *IT
            RES(J,L) = S(J)*BPLUS(L) + RFLUX(J,L)*EPLUS(L)
            DO 200 K = 1, NPDE
               RES(J,L-1) = RES(J,L-1) - (APLUS(L)*PF(J,K)+PL(J,K))
     *                      /(APLUS(L)+BPLUS(L-1))*UDOT(K,L-1)
               PL(J,K) = BPLUS(L)*PF(J,K)
  200       CONTINUE
            CS(J) = CF(J)
            DS(J) = DC(J)
  220    CONTINUE
  240 CONTINUE
C
C ... Finish processing boundary conditions at X(1) ...
C
      DO 260 J = 1, NPDE
         RES(J,1) = G(J)*IT
  260 CONTINUE
C
C ... Finally, process the right-hand boundary conditions ...
C
      IBND = 1
      CALL BNDR(PLFBND,BNDPFF,T,G,U,UB,X,NPDE,NPTS,IBND,NV,V,VDOT,IZ(3),
     *          IIFLAG)
C
      IF (IIFLAG.EQ.2) THEN
         IRES = IZ(3)
         RETURN
      END IF
C
      IF (NVIT.EQ.0) THEN
         CALL BNDR(PLFBND,BNDPFF,T,BT,U,UB,X,NPDE,NPTS,IBND,NV,V,VT,
     *             IZ(3),IIFLAG)
C
         IF (IIFLAG.EQ.2) THEN
            IRES = IZ(3)
            RETURN
         END IF
C
         DO 280 J = 1, NPDE
            G(J) = G(J) - BT(J,1)
  280    CONTINUE
      END IF
      DO 300 J = 1, NPDE
         RES(J,NPTS) = G(J)*IT
  300 CONTINUE
C
C ... Illegal problem spec encountered ...
C
      IF (ABS(TEMP).GE.UROUND) THEN
         IIFLAG = 1
         RETURN
      END IF
      DO 320 I = 1, 3
         IF (IZ(I).NE.1) THEN
            IR = IZ(I)
            GO TO 340
         END IF
  320 CONTINUE
  340 CONTINUE
      RETURN
      END
