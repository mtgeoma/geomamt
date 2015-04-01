      SUBROUTINE D03PKR(NEQN,NPDE,NPTS,X,U,WK,IWK,M,TS,IBAND,ITIME,
     *                  REMESH,MAXNPT,NV,NXI,XI,XFIX,NXFIX,IXFIX,NIXFIX,
     *                  UVINIT,PKFPDE,PDEPEF,SPDEF1,MONFFD,MONFKB)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C   --------------------------------------------------------------------
C
C     K.B. VERSION OF D03PHR
C
C   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      In the case when mixed ODE/PDE problems are
C      solved and the body of the PDE equation (one of
C      the functions C,F or R in the problem spec above )
C      depends on the ODE variables , or the ODE
C      variables depend on non boundary PDE values the
C      normal banded structure of the ODE system is
C      destroyed . When solving mixed ODE /PDE
C      problems sparse or full linear algebra routines
C      should be used in SPRINT.
C
C     IBAND ; is set to NEQN if NV > 0.
C
C     ITIME ;  should be set to one on the first call . On this call
C            the user supplied mesh X(NPTS) is put in the first
C            NPTS locations of work. On a subsequent call the mesh
C            used , which may have been changed by a mesh
C            modification monitor routine, may be put into X(NPTS)
C            calling INITSK with ITIME = 2.
C            This parameter is set to -1 if an error is found by
C            this routine.
C
C     SEVEN IMPORTANT PARAMETERS ARE PASSED ACROSS FROM HERE IN
C     COMMON /PD03PC/  NNPDE, NNPTS, MM, NNV, NNXI, NVST, NPMAX
C
C     Detailed description of workspace :
C
C     The workspace WK(IWK) is used to pass arrays and vectors to
C     the routines  D03PKJ , D03PZW and D03PEN .
C
C     size  :  IWK must be >=   NPDE*(MAXNPT + NXI*6 + 15 + 2*NPDE)
C                                + 7*MAXNPT + 2*NXFIX  + NXI + NV
C
C     Structure        Name in code   Purpose in code is to hold
C
C     WK(1) - WK(I1-1)            array of size (13,NPDE) used by the
C                                 D03PHL discretisation routine
C
C       ... The following parts of the workspace are used ...
C       ... in semi-discretising mixed ODE/PDE problems   ...
C
C     WK(I1) - WK(I2-1) UI       array UI(NPDE,NXI) used to hold the
C                                PDE components at the coupling pts
C     WK(I2) - WK(I3-1) UXI      space derivs corressponding to UI
C     WK(I3) - WK(I4-1) RI       flux corress to UI array
C     WK(I4) - WK(I5-1) UTI      time deriv corressponding to UI
C     WK(I5) - WK(I6-1) UTXI     space deriv of array UTI
C     WK(I6) - WK(I7-1) XI       Coupling points to link PDE to ODE
C     WK(I7) - WK(I8-1) VDUM     array of zeroes
C
C      ... The following three parts of WK hold the mesh ...
C      ... point information used in semi-discretisation ...
C
C     WK(I8) - WK(I9-1)  X(NPTS)   spatial mesh points used by the code
C
C     WK(I9) - WK(I10-1) R(NPDE,NPTS+1) flux values at the boundaries
C                         and at the points midway between the mesh
C                         points. used in mesh modification only
C
C     WK(I10) - WK(I14)   arrays used in the remeshing process to
C                         generate the new mesh
C
C     WK(I14) - WK(I16)   Arrays used in applying remeshing to problems
C                         with material interfaces.
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TS
      INTEGER           IBAND, ITIME, IWK, M, MAXNPT, NEQN, NIXFIX,
     *                  NPDE, NPTS, NV, NXFIX, NXI
      LOGICAL           REMESH
C     .. Array Arguments ..
      DOUBLE PRECISION  U(*), WK(IWK), X(NPTS), XFIX(*), XI(*)
      INTEGER           IXFIX(*)
C     .. Subroutine Arguments ..
      EXTERNAL          MONFFD, MONFKB, PDEPEF, PKFPDE, SPDEF1, UVINIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  CONST, STPRAT, TDUM
      INTEGER           I1, I10, I11, I12, I13, I14, I15, I16, I2, I3,
     *                  I4, I5, I6, I7, I8, I9, IIFLAG, J1, J10, J11,
     *                  J2, J3, J4, J5, J6, J7, J8, J9, JPMINF, K2, K3,
     *                  MM, NNPDE, NNPTS, NNV, NNXFIX, NNXI, NOPUSR,
     *                  NPMAX, NVST
      LOGICAL           REMSH1
      CHARACTER*6       PDCODE, RMTYPE
C     .. Arrays in Common ..
      INTEGER           IDUM(2)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, D03PKH
C     .. Common blocks ..
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /GD03PC/CONST, STPRAT, JPMINF, NOPUSR
      COMMON            /HD03PC/TDUM, REMSH1, IDUM
      COMMON            /PD03PC/NNPDE, NNPTS, MM, NNV, NNXI, NVST,
     *                  NPMAX, NNXFIX
      COMMON            /QD03PC/I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I11, I12, I13, I14, I15, I16
      COMMON            /RD03PC/K2, K3
      COMMON            /SD03PC/J1, J2, J3, J4, J5, J6, J7, J8, J9, J10,
     *                  J11
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /PD03PC/, /FD03PC/, /QD03PC/, /GD03PC/,
     *                  /SD03PC/, /HD03PC/, /RD03PC/, /XD03PC/
C     .. Executable Statements ..
C
      IF (ITIME.GT.1) THEN
C
C        ... Update the number and positioning of the mesh points ...
C
         NPTS = NNPTS
         DO 20 I = 1, NPTS
            X(I) = WK(I8+I-1)
   20    CONTINUE
         RETURN
      ELSE
C
C        ... Zero the workspace to prevent undefined ...
C        ... variables being used.                   ...
C
         DO 40 I = 1, IWK
            WK(I) = 0.0D0
   40    CONTINUE
      END IF
      MM = M
      NNPDE = NPDE
C
C     ... If ITIME>1, NPTS and NEQN updated in D03PEN ...
C     ... remeshing monitor subroutine.               ...
C
      NNPTS = NPTS
      NEQN = NPDE*NPTS + NV
      NPMAX = MAXNPT
      NNXFIX = NXFIX
C
      J2 = 1 + NPDE
      J3 = J2 + NPDE
      J4 = J3 + NPDE
      J5 = J4 + NPDE
      J6 = J5 + NPDE
      J7 = J6 + NPDE
      J8 = J7 + NPDE*NPDE
      J9 = J8 + NPDE*NPDE
      J10 = J9 + NPDE
      J11 = J10 + NPDE*NPDE
      I1 = (14+3*NPDE)*NPDE + 1
      I2 = I1 + NPDE*NXI
      I3 = I2 + NPDE*NXI
      I4 = I3 + NPDE*NXI
      I5 = I4 + NPDE*NXI
      I6 = I5 + NPDE*NXI*2
      I7 = I6 + NXI
      I8 = I7 + NV
      I9 = I8 + MAXNPT
      I10 = I9 + NPDE*(MAXNPT+1)
      I11 = I10 + MAXNPT
      I12 = I11 + MAXNPT
      I13 = I12 + MAXNPT
      K2 = I13 + MAXNPT
      K3 = I13 + MAXNPT*2
C
      IF (NXFIX.EQ.0) THEN
         I14 = IWK
         I15 = IWK
      ELSE
         I14 = I13 + MAXNPT*3
         I15 = I14 + NXFIX
      END IF
C
      I = I8 + NPDE + (NPDE+7)*MAXNPT + 2*NXFIX - 1
      REMSH1 = REMESH
      NNV = NV
      NNXI = NXI
      NVST = NPDE*NPTS + 1
      IF (NV.EQ.0) NVST = NVST - 1
C
      DO 60 I = 1, NXI
         WK(I6+I-1) = XI(I)
   60 CONTINUE
C
C     ... Initialise ODE and PDE variables by appropriate calls ...
C
      CALL UVINIT(NPDE,NPTS,NXI,X,XI,U,NV,U(NVST))
C
      DO 80 I = 1, NPTS
         WK(I8+I-1) = X(I)
   80 CONTINUE
C
C     ... Call to the routine to adjust the initial mesh ...
C
      NNXFIX = NXFIX
      IF (REMESH) THEN
         DO 100 I = 1, NXFIX
            WK(I14+I-1) = XFIX(I)
  100    CONTINUE
C
         I = 1
         J = 1
C
         CALL D03PKH(PKFPDE,PDEPEF,TS,M,U,NPDE,NNPTS,WK,IWK,U(NVST),NV,
     *               NPMAX,I,XI,NXI,WK(I14),IXFIX,NNXFIX,J,UVINIT,
     *               SPDEF1,MONFFD,MONFKB)
C
         IF (IIFLAG.EQ.1 .OR. IIFLAG.EQ.2) THEN
            CALL D02NNQ(
     *' User-supplied routine PDEDEF sets illegal value of
     * IRES (=I1) during initial remeshing.',1,1,I,0,0,0.0D0,0.0D0)
            RETURN
         END IF
C
         IF (I.NE.1 .OR. J.LT.0) THEN
            IF (I.NE.1) THEN
               CALL D02NNQ(
     *' User-supplied routine PDEDEF sets IRES (=I1) to 2 or 3
     * during initial remeshing - check problem specification.',1,1,I,0,
     *                     0,0.0D0,0.0D0)
               IIFLAG = 2
            ELSE IF (J.LT.0) THEN
C              CALL D02NNQ(
C              *' Remesh routines indicated that calculated or supplied
C              * spatial mesh did not match discontinuity points
C              * supplied in the array XFIX ',1,0,0,0,0,0.0D0,0.0D0)
               ITIME = -1
               RETURN
            END IF
         END IF
C
C        ... Resets NNPTS and NEQN after call to D03PHH since ...
C        ... this may change NPTS.                            ...
C
         NEQN = NPDE*NNPTS + NV
         NPTS = NNPTS
         NVST = NPDE*NPTS
         IF (NV.GT.0) NVST = NPDE*NPTS + 1
C
         CALL UVINIT(NPDE,NPTS,NXI,WK(I8),XI,U,NV,U(NVST))
C
C        ... Update the mesh returned to the user ...
C
         DO 120 I = 1, NPTS
            X(I) = WK(I8+I-1)
  120    CONTINUE
      END IF
C
      PDCODE = 'KELLER'
      RETURN
      END
