      SUBROUTINE D03PDQ(PDEFCN,BNDARY,DUMPD1,DUMPD2,NEQN,T,U,UDOT,RESD,
     *                  IRES,WK,IWK,SPDEFN,SBNDR,SODEFN)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ----------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This is the chebyshev global element routine to evaluate the
C     residual of the implicit set of ODEs defined by
C
C        RESIDUAL  =  A(u,t)*du/dt  -  F(u,t)
C
C                     PARAMETER LIST
C                    ----------------
C
C     NEQN(1) = N     Number of ODEs in time.
C     M               Polar parameter for PDEs.
C     NPDE            Number of parabolic equations.
C     NPTS            Number of spatial grid points.
C     T               Current time integration level, > 0.0.
C     U(N)            Current solution vector.
C     RESD(N)         Vector which will contain the residual on exit
C     UDOT(N)         Current estimate of du/dt
C     WK(IWK)         Workspace - defined in inital
C     IRES            Indicator for residual routine.
C                     On entry = -1 then evaluate those parts of the
C                                   residual only which involve du/dt
C                              =  0 Evaluate the full residual.
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ----------------------------------------------------------------------
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, IWK, NEQN
C     .. Array Arguments ..
      DOUBLE PRECISION  RESD(*), U(*), UDOT(*), WK(IWK)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, DUMPD1, DUMPD2, PDEFCN, SBNDR, SODEFN,
     *                  SPDEFN
C     .. Scalars in Common ..
      INTEGER           I10, I10A, I10B, I11, I11A, I11B, I12, I13, I14,
     *                  I15, I16, I17, I18, I19, I2, I3, I4, I5, I6, I7,
     *                  I8, I9, IDEV, IIFLAG, ITRACE, M, NEL, NPDE,
     *                  NPTL, NPTS, NV, NVST, NXI
C     .. Local Scalars ..
      INTEGER           I, IBK, IFL, IJ, IR, ITYPE, IV, J, K, N
C     .. External Subroutines ..
      EXTERNAL          D03PDR, D03PDS
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /ED03PD/NEL, NPTL, NPDE, NPTS, M, NV, NXI, NVST
      COMMON            /FD03PD/I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I10A, I10B, I11, I11A, I11B, I12, I13, I14, I15,
     *                  I16, I17, I18, I19
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /ED03PD/, /FD03PD/, /AD02NM/
C     .. Executable Statements ..
C
      N = NEQN
      DO 20 J = 1, N
         RESD(J) = 0.0D0
   20 CONTINUE
      IIFLAG = 0
      IBK = NEL + 1
      IV = NPTS*NPDE
      IF (NV.GT.0) THEN
         IV = NVST
         IFL = 0
         IF (NXI.GT.0) THEN
C
C ... Generate the solution values space derivs and fluxes at the ...
C ... coupling points.                                            ...
C
            ITYPE = 3
C
            CALL D03PDR(PDEFCN,DUMPD1,NXI,WK(I18),WK(I13),ITYPE,U,NPTS,
     *                  NPDE,NEL,NPTL,WK,WK(I10),WK(I5),IBK,IFL,NV,U(IV)
     *                  ,UDOT(IV),WK(I11),T,IR,SPDEFN)
C
            IF (IIFLAG.EQ.2) THEN
               IRES = IR
               GO TO 160
            END IF
C
            IF (IR.NE.1) GO TO 140
C
C ... Generate those parts of fluxes that depend on time ...
C ... derivitives.                                       ...
C
            IF (IRES.EQ.-1) THEN
               ITYPE = 4
               DO 60 J = 1, NXI
                  DO 40 K = 1, NPDE
                     IJ = NPDE*(J-1) + K - 1
                     WK(I16+IJ) = WK(I15+IJ)
   40             CONTINUE
   60          CONTINUE
C
               CALL D03PDR(PDEFCN,DUMPD1,NXI,WK(I18),WK(I13),ITYPE,U,
     *                     NPTS,NPDE,NEL,NPTL,WK,WK(I10),WK(I5),IBK,IFL,
     *                     NV,U(IV),WK(I19),WK(I11),T,IR,SPDEFN)
C
               IF (IIFLAG.EQ.2) THEN
                  IRES = IR
                  GO TO 160
               END IF
C
               IF (IR.NE.1) GO TO 140
C
               DO 100 J = 1, NXI
                  DO 80 K = 1, NPDE
                     IJ = NPDE*(J-1) + K - 1
                     WK(I15+IJ) = WK(I16+IJ) - WK(I15+IJ)
   80             CONTINUE
  100          CONTINUE
            END IF
C
C ... Generate time deriv values and their space derivs at the ...
C ... coupling points.                                         ...
C
            ITYPE = 2
C
            CALL D03PDR(PDEFCN,DUMPD1,NXI,WK(I18),WK(I16),ITYPE,UDOT,
     *                  NPTS,NPDE,NEL,NPTL,WK,WK(I10),WK(I5),IBK,IFL,NV,
     *                  U(IV),UDOT(IV),WK(I11),T,IR,SPDEFN)
C
            IF (IIFLAG.EQ.2) THEN
               IRES = IR
               GO TO 160
            END IF
C
            IF (IR.NE.1) GO TO 140
         END IF
C
C ... Call the routine to define the auxillary ode residual ...
C
         CALL SODEFN(NPDE,T,NV,U(IV),UDOT(IV),NXI,WK(I18),WK(I13),
     *               WK(I14),WK(I15),WK(I16),WK(I17),RESD(IV),IRES)
C
         IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
            IIFLAG = 2
            RETURN
         END IF
C
         IF (ABS(IRES).NE.1) THEN
            RETURN
         END IF
      END IF
C
C ... Call the co collocation discretisation routine ...
C
      IR = 1
C
      CALL D03PDS(PDEFCN,BNDARY,DUMPD1,DUMPD2,NPDE,NPTS,T,U,RESD,UDOT,M,
     *            WK(I6),WK,WK(I2),WK(I5),WK(I7),WK(I8),WK(I9),WK(I10),
     *            WK(I11),NEL,NPTL,WK(I4),WK(I12),IRES,WK(I10A),WK(I11A)
     *            ,WK(I11B),WK(I10B),NV,U(IV),UDOT(IV),WK(I19),SPDEFN,
     *            SBNDR)
C
      IF (IIFLAG.EQ.1 .OR. IIFLAG.EQ.2) GO TO 160
C
      DO 120 I = 1, N
         RESD(I) = -RESD(I)
  120 CONTINUE
      RETURN
C
  140 CONTINUE
      IRES = IR
  160 CONTINUE
      RETURN
      END
