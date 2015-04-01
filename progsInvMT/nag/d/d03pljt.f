      SUBROUTINE D03PLJ(PLFPDE,PLFBND,PDEPFF,BNDPFF,NEQN,T,U,UDOT,FF,
     *                  IRES,RWK,NRWK,PDEFN,FLXPFF,FLXPLF,BNDY,ODEFN)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C----------------------------------------------------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ODE residual routine that performs spatial discretisation
C  of mixed ODE/PDE problems in one space dimension.
C  The method used is the Skeel finite difference discretisation
C  with upwinding of convective terms.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, NEQN, NRWK
C     .. Array Arguments ..
      DOUBLE PRECISION  FF(*), RWK(NRWK), U(*), UDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPFF, BNDY, FLXPFF, FLXPLF, ODEFN, PDEFN,
     *                  PDEPFF, PLFBND, PLFPDE
C     .. Scalars in Common ..
      INTEGER           I1, I10, I11, I12, I13, I14, I15, I16, I2, I3,
     *                  I4, I5, I6, I7, I8, I9, ID, IIFLAG, IRNITE, J1,
     *                  J10, J11, J12, J13, J14, J15, J16, J17, J18, J2,
     *                  J3, J4, J5, J6, J7, J8, J9, K2, K3, M, NPDE,
     *                  NPTS, NV, NVST, NXFIX, NXI
      CHARACTER*6       PDCODE, RMTYPE
C     .. Local Scalars ..
      INTEGER           I, IFL, IR, IT, IV, NRST, NULST, NURST
C     .. External Subroutines ..
      EXTERNAL          D03PLK, D03PLL, D03PZW
C     .. Common blocks ..
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /PD03PC/NPDE, NPTS, M, NV, NXI, NVST, ID, NXFIX
      COMMON            /QD03PC/I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I11, I12, I13, I14, I15, I16
      COMMON            /RD03PC/K2, K3
      COMMON            /UD03PC/IRNITE
      COMMON            /XD03PC/IIFLAG
      COMMON            /ZD03PC/J1, J2, J3, J4, J5, J6, J7, J8, J9, J10,
     *                  J11, J12, J13, J14, J15, J16, J17, J18
C     .. Save statement ..
      SAVE              /PD03PC/, /FD03PC/, /QD03PC/, /ZD03PC/,
     *                  /UD03PC/, /RD03PC/, /XD03PC/
C     .. Executable Statements ..
C
      IF (NV.GT.0) THEN
         DO 20 I = 1, NV
            RWK(I7+I-1) = 0.0D0
   20    CONTINUE
      END IF
      IR = 1
C
C  call the D03PLK convective flux routine..
C
C  calculate start of convective flux array for D03PLK ..
      NRST = NRWK - 2*NPDE - NPDE*NPTS + 1
C  calculate start of ULEFT and URIGHT arrays for D03PLK ..
      NULST = NRST + NPDE*NPTS
      NURST = NULST + NPDE
C
      CALL D03PLK(NPDE,NPTS,RWK(I8),T,U,RWK(J3),RWK(NULST),RWK(NURST),
     *            NV,U(NVST),FLXPFF,FLXPLF,RWK(NRST),IR,IRES)
C
      IF (IIFLAG.EQ.1 .OR. IIFLAG.EQ.2) THEN
         RETURN
      END IF
      IF (IR.NE.1) THEN
         IRES = IR
         RETURN
      END IF
C
C  call the D03PLL discretisation routine using partitioned workspace ..
C
      CALL D03PLL(PLFPDE,PLFBND,PDEPFF,BNDPFF,RWK(I8),T,U,UDOT,NPDE,
     *            NPTS,FF,RWK,RWK(J2),RWK(J3),RWK(J4),RWK(J5),RWK(J6),
     *            RWK(J7),RWK(J8),IRES,U(NVST),UDOT(NVST),NV,RWK(J9),
     *            RWK(J10),RWK(J11),RWK(J12),RWK(J13),RWK(J14),RWK(J15),
     *            RWK(J16),RWK(J17),RWK(J18),RWK(I7),IR,RWK(NRST),PDEFN,
     *            BNDY,RWK(I9),RWK(I10),RWK(I11),RWK(I12),RWK(I13),
     *            RWK(K2),RWK(K3))
C
      IF (IIFLAG.EQ.1 .OR. IIFLAG.EQ.2) THEN
         RETURN
      END IF
      IF (IR.NE.1) THEN
         IRES = IR
         RETURN
      END IF
      IFL = 0
      IV = NVST
      IF (NV.GT.0) THEN
         IF (NXI.GT.0) THEN
C
C ... Generate U values at coupling points ...
C
            IT = 2
            CALL D03PZW(RWK(I6),RWK(I1),NXI,RWK(I8),M,U,NPTS,NPDE,IT,
     *                  IFL)
C
C ... Generate time derivatives at coupling points ...
C
            IT = 1
            CALL D03PZW(RWK(I6),RWK(I4),NXI,RWK(I8),M,UDOT,NPTS,NPDE,IT,
     *                  IFL)
C
         END IF
C
C ... Define the auxillary ODE residual. ...
C
         CALL ODEFN(NPDE,T,NV,U(IV),UDOT(IV),NXI,RWK(I6),RWK(I1),RWK(I2)
     *              ,RWK(I4),FF(IV),IRES)
C
         IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
            IIFLAG = 2
            RETURN
         END IF
C
      END IF
      RETURN
      END
