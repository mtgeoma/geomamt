      SUBROUTINE D03PKJ(PKFPDE,PKFBND,PDEPEF,BNDPEF,NEQN,T,U,UDOT,RES,
     *                  IRES,WK,NWKRES,PDEFN,BNDY,ODEFN)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C-----------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     D03P   standard discretisation routine for Keller Box method
C     developed from that of Ron Furzeland, 1987.
C     unknowns are U(NPDE,NPTS) held in 1-D array U.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, NEQN, NWKRES
C     .. Array Arguments ..
      DOUBLE PRECISION  RES(*), U(*), UDOT(*), WK(NWKRES)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPEF, BNDY, ODEFN, PDEFN, PDEPEF, PKFBND,
     *                  PKFPDE
C     .. Scalars in Common ..
      INTEGER           I1, I10, I11, I12, I13, I14, I15, I16, I2, I3,
     *                  I4, I5, I6, I7, I8, I9, ID, IIFLAG, IONE,
     *                  IRNITE, IW1, IW2, IW3, IW4, IW5, IW6, IW7, IW8,
     *                  J1, J10, J11, J2, J3, J4, J5, J6, J7, J8, J9,
     *                  K2, K3, M, MM, NLEFT, NPDE, NPTS, NRIGHT, NV,
     *                  NVST, NXFIX, NXI
      CHARACTER*6       PDCODE, RMTYPE
C     .. Local Scalars ..
      DOUBLE PRECISION  HX, XMEAN
      INTEGER           I, IB, IBNDY, IFL, IJ, IJP1, IJST, IRET, IT, J,
     *                  NPM
C     .. Local Arrays ..
      DOUBLE PRECISION  DV(1)
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, D03PZW
C     .. Common blocks ..
      COMMON            /BD03PC/NLEFT, NRIGHT, MM
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /PD03PC/NPDE, NPTS, M, NV, NXI, NVST, ID, NXFIX
      COMMON            /QD03PC/I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I11, I12, I13, I14, I15, I16
      COMMON            /RD03PC/K2, K3
      COMMON            /SD03PC/J1, J2, J3, J4, J5, J6, J7, J8, J9, J10,
     *                  J11
      COMMON            /TD03PC/IONE, IW1, IW2, IW3, IW4, IW5, IW6, IW7,
     *                  IW8
      COMMON            /UD03PC/IRNITE
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /BD03PC/, /FD03PC/, /PD03PC/, /QD03PC/,
     *                  /RD03PC/, /SD03PC/, /TD03PC/, /UD03PC/, /XD03PC/
C     .. Executable Statements ..
C
      NPM = NPTS - 1
C
      IONE = I8
C                IONE is start address of the MESH
      IW1 = 1
C                IW1  is start address of UMEAN array
      IW2 = J2
C                IW2  is start address of UDOT mean array
      IW3 = J3
C                IW3  is start address of DUDX mean array
      IW4 = I6
C                IW4  is start address of XI array
      IW5 = I1
C                IW5  is start address of UI array
      IW6 = I2
C                IW6  is start address of UXI array
      IW7 = I3
C                IW7  is start address of UTXI array
      IW8 = I4
C
C     ... Just left hand side required for IRES=-1 ...
C     ... Sets up residuals        RES(1),.....RES(NLEFT) ...
      IF (NLEFT.EQ.0) GO TO 20
C
      IBNDY = 0
C
      IF (NV.EQ.0) THEN
         DV(1) = 0.0D0
C      CALL BNDY(PKFBND,BNDPEF,T,WK(IONE),IBNDY,NPDE,U,UDOT,NV,DV,
C     *          DV,NLEFT,RES,IRES,IIFLAG)
         CALL BNDY(PKFBND,BNDPEF,T,IBNDY,NPDE,U,UDOT,NV,DV,DV,NLEFT,RES,
     *             IRES,IIFLAG)
      ELSE
C      CALL BNDY(PKFBND,BNDPEF,T,WK(IONE),IBNDY,NPDE,U,UDOT,NV,U(NVST),
C     *          UDOT(NVST),NLEFT,RES,IRES,IIFLAG)
         CALL BNDY(PKFBND,BNDPEF,T,IBNDY,NPDE,U,UDOT,NV,U(NVST),
     *             UDOT(NVST),NLEFT,RES,IRES,IIFLAG)
      END IF
C
C     ... Coupled ODES (SUBR.ODEFN) are held in ...
C     ... U(NVST=NPDE*NPTS+1) ...
C
      IF (IRES.EQ.2 .OR. IRES.EQ.3) THEN
         CALL D02NNQ(
     *' Error in residual routine D03PKJ - user routine BNDARY
     *  has signalled non-zero error return code (=I1) ',1,0,IRES,0,0,
     *               0.0D0,0.0D0)
         RETURN
      END IF
C
   20 IB = NLEFT + 1
C
      DO 60 I = 1, NPM
         HX = WK(I+IONE) - WK(I+IONE-1)
         XMEAN = (WK(I+IONE)+WK(I+IONE-1))/2.0D0
         IJST = (I-1)*NPDE
C
         DO 40 J = 1, NPDE
            IJ = IJST + J
            IJP1 = IJ + NPDE
            WK(IW1+J-1) = (U(IJP1)+U(IJ))/2.0D0
            WK(IW2+J-1) = (UDOT(IJP1)+UDOT(IJ))/2.0D0
            WK(IW3+J-1) = (U(IJP1)-U(IJ))/HX
   40    CONTINUE
C
         IF (NV.EQ.0) THEN
            DV(1) = 0.0D0
            CALL PDEFN(PKFPDE,PDEPEF,T,XMEAN,NPDE,WK(IW1),WK(IW2),
     *                 WK(IW3),NV,DV,DV,RES(IB),IRES,IIFLAG)
         ELSE
            CALL PDEFN(PKFPDE,PDEPEF,T,XMEAN,NPDE,WK(IW1),WK(IW2),
     *                 WK(IW3),NV,U(NVST),UDOT(NVST),RES(IB),IRES,
     *                 IIFLAG)
         END IF
C
         IF (IRES.EQ.2 .OR. IRES.EQ.3) THEN
C
            CALL D02NNQ(
     *' Error in residual routine D03PKJ - user routine PDEDEF
     *  has signalled non-zero error return code (=I1) ',1,0,IRET,0,0,
     *                  0.0D0,0.0D0)
C
            RETURN
         END IF
C
         IB = IB + NPDE
C
   60 CONTINUE
C
C     ... Right-hand boundary   (NRIGHT=NPDE-NLEFT) ...
C     ... Sets up residuals     RES(NEQ-NRIGHT+1), ..., RES(NEQ) ...
C
      IF (NRIGHT.EQ.0) GO TO 80
      IJST = (NPTS-1)*NPDE + 1
      IBNDY = 1
      I = IONE + NPTS - 1
C
      IF (NV.EQ.0) THEN
         DV(1) = 0.0D0
C      CALL BNDY(PKFBND,BNDPEF,T,WK(I),IBNDY,NPDE,U(IJST),UDOT(IJST),NV,
C     *          DV,DV,NRIGHT,RES(IB),IRES,IIFLAG)
         CALL BNDY(PKFBND,BNDPEF,T,IBNDY,NPDE,U(IJST),UDOT(IJST),NV,DV,
     *             DV,NRIGHT,RES(IB),IRES,IIFLAG)
      ELSE
C      CALL BNDY(PKFBND,BNDPEF,T,WK(I),IBNDY,NPDE,U(IJST),UDOT(IJST),NV,
C     *          U(NVST),UDOT(NVST),NRIGHT,RES(IB),IRES,IIFLAG)
         CALL BNDY(PKFBND,BNDPEF,T,IBNDY,NPDE,U(IJST),UDOT(IJST),NV,
     *             U(NVST),UDOT(NVST),NRIGHT,RES(IB),IRES,IIFLAG)
      END IF
C
      IF (IRES.EQ.2 .OR. IRES.EQ.3) THEN
         CALL D02NNQ(
     *' Error in residual routine D03PKJ - user routine BNDARY
     *  has signalled non-zero error return code (=I1) ',1,0,IRES,0,0,
     *               0.0D0,0.0D0)
         RETURN
      END IF
C
   80 IFL = 0
C
      IF (NV.GT.0) THEN
         IF (NXI.GT.0) THEN
            IT = 2
C
            CALL D03PZW(WK(IW4),WK(IW5),NXI,WK(IONE),M,U,NPTS,NPDE,IT,
     *                  IFL)
C
            CALL D03PZW(WK(IW4),WK(IW7),NXI,WK(IONE),M,UDOT,NPTS,NPDE,
     *                  IT,IFL)
         END IF
C
         CALL ODEFN(NPDE,T,NV,U(NVST),UDOT(NVST),NXI,WK(IW4),WK(IW5),
     *              WK(IW6),WK(IW7),RES(NVST),IRES)
C
      END IF
C
      RETURN
      END
