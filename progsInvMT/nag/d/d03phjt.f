      SUBROUTINE D03PHJ(PHFPDE,PHFBND,PDEPCF,BNDPCF,NEQN,T,U,UDOT,FF,
     *                  IRES,RWK,NRWK,PDEFN,BNDY,ODEFN)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C----------------------------------------------------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  ODE residual routine that performs spatial discretisation
C  of mixed ODE / PDE. problems in one space dimension.
C  the method used is the Skeel  finite difference discretisation.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, NEQN, NRWK
C     .. Array Arguments ..
      DOUBLE PRECISION  FF(*), RWK(NRWK), U(*), UDOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPCF, BNDY, ODEFN, PDEFN, PDEPCF, PHFBND,
     *                  PHFPDE
C     .. Scalars in Common ..
      INTEGER           I1, I10, I11, I12, I13, I14, I15, I16, I2, I3,
     *                  I4, I5, I6, I7, I8, I9, ID, IIFLAG, IRNITE, J1,
     *                  J10, J11, J2, J3, J4, J5, J6, J7, J8, J9, K2,
     *                  K3, M, NPDE, NPTS, NV, NVST, NXFIX, NXI
      CHARACTER*6       PDCODE, RMTYPE
C     .. Local Scalars ..
      INTEGER           I, IFL, IR, IT, IV, N
C     .. External Subroutines ..
      EXTERNAL          D03PHK, D03PHL, D03PZW
C     .. Common blocks ..
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /PD03PC/NPDE, NPTS, M, NV, NXI, NVST, ID, NXFIX
      COMMON            /QD03PC/I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I11, I12, I13, I14, I15, I16
      COMMON            /RD03PC/K2, K3
      COMMON            /SD03PC/J1, J2, J3, J4, J5, J6, J7, J8, J9, J10,
     *                  J11
      COMMON            /UD03PC/IRNITE
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /PD03PC/, /FD03PC/, /QD03PC/, /SD03PC/,
     *                  /UD03PC/, /RD03PC/
C     .. Executable Statements ..
C
C
C  CALL THE D03PHL DISCRETISATION ROUTINE USING PARTITIONED WORKSPACE.
C
      IF (NV.GT.0) THEN
         DO 20 I = 1, NV
            RWK(I7+I-1) = 0.0D0
   20    CONTINUE
      END IF
      IR = 1
C
      CALL D03PHL(PHFPDE,PHFBND,PDEPCF,BNDPCF,RWK(I8),T,U,UDOT,NPDE,
     *            NPTS,M,FF,RWK,RWK(J2),RWK(J3),RWK(J4),RWK(J5),RWK(J6),
     *            RWK(J7),RWK(J8),RWK(I9),IRES,U(NVST),UDOT(NVST),NV,
     *            RWK(J9),RWK(J10),RWK(J11),RWK(I7),IR,PDEFN,BNDY,
     *            RWK(I9),RWK(I10),RWK(I11),RWK(I12),RWK(I13),RWK(K2),
     *            RWK(K3))
C
      IF (IIFLAG.EQ.1 .OR. IIFLAG.EQ.2) THEN
         RETURN
      END IF
      IF (IR.NE.1) THEN
         IRES = IR
         RETURN
      END IF
      N = NEQN
      IFL = 0
      IV = NVST
      IF (NV.GT.0) THEN
         IF (NXI.GT.0) THEN
            IT = 2
C
C ... Generate U  values and space derivs at coupling points ...
C
            CALL D03PZW(RWK(I6),RWK(I1),NXI,RWK(I8),M,U,NPTS,NPDE,IT,
     *                  IFL)
C
C ... Generate flux values at coupling points ...
C
            CALL D03PHK(RWK(I6),RWK(I3),NXI,RWK(I8),RWK(I9),NPTS,NPDE,
     *                  IFL)
C
C ... Generate time derivs and their space derivs at coupling pts ...
C
            CALL D03PZW(RWK(I6),RWK(I4),NXI,RWK(I8),M,UDOT,NPTS,NPDE,IT,
     *                  IFL)
C
         END IF
C
C ... Define the auxillary ODE residual. ...
C
         CALL ODEFN(NPDE,T,NV,U(IV),UDOT(IV),NXI,RWK(I6),RWK(I1),RWK(I2)
     *              ,RWK(I3),RWK(I4),RWK(I5),FF(IV),IRES)
C
         IF ((IRES.LT.-1) .OR. (IRES.EQ.0) .OR. (IRES.GT.3)) THEN
            IIFLAG = 2
            RETURN
         END IF
C
      END IF
      RETURN
      END
