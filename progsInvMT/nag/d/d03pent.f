      SUBROUTINE D03PEN(PKFPDE,PKFBND,PDEPEF,BNDPEF,NEQN,T,HLAST,H,Y,
     *                  YDOT,YSAVE,NYH,R,ACOR,RESWK,NRESWK,WKMON,NWKMON,
     *                  IMON,INLN,HMIN,HMAX,MONFFD,MONFKB,PDEFN,BNDR,
     *                  IRES,IXFIX,NIXFIX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C ----------------------------------------------------------------------
C     MONITOR routine for remeshing (Keller Box Scheme)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HLAST, HMAX, HMIN, T
      INTEGER           IMON, INLN, IRES, NEQN, NIXFIX, NRESWK, NWKMON,
     *                  NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(*), R(*), RESWK(NRESWK), WKMON(NWKMON),
     *                  Y(*), YDOT(*), YSAVE(NYH,*)
      INTEGER           IXFIX(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPEF, BNDR, MONFFD, MONFKB, PDEFN, PDEPEF,
     *                  PKFBND, PKFPDE
C     .. Scalars in Common ..
      DOUBLE PRECISION  CONST, DUNFLO, DXMESH, HOLD, STPRAT, TOUT,
     *                  UROUND
      INTEGER           I1, I10, I11, I12, I13, I14, I15, I16, I2, I3,
     *                  I4, I5, I6, I7, I8, I9, ICOUNT, IDEV, IOVFLO,
     *                  IPMINF, ITRACE, J1, J10, J11, J2, J3, J4, J5,
     *                  J6, J7, J8, J9, K2, K3, KCUR, M, MAXNPT, NINTER,
     *                  NOPUSR, NPDE, NPTS, NRMESH, NV, NVST, NXFIX, NXI
      LOGICAL           REMESH
      CHARACTER*6       PDCODE, RMTYPE
C     .. Arrays in Common ..
      INTEGER           NDUM(6)
C     .. Local Scalars ..
      INTEGER           I, IFAIL, K, MIDNPT, N, NF, NIP, NOLD, NOP
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          D02NNN, D03PRH, D03PRK, D03PRM, D03PRP, D03PRQ,
     *                  X04ABF, X04BAF
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /BD02NM/HOLD, NDUM, NINTER, KCUR
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /GD03PC/CONST, STPRAT, IPMINF, NOPUSR
      COMMON            /HD03PC/TOUT, REMESH, NRMESH, ICOUNT
      COMMON            /MD03PC/DXMESH
      COMMON            /PD03PC/NPDE, NPTS, M, NV, NXI, NVST, MAXNPT,
     *                  NXFIX
      COMMON            /QD03PC/I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I11, I12, I13, I14, I15, I16
      COMMON            /RD03PC/K2, K3
      COMMON            /SD03PC/J1, J2, J3, J4, J5, J6, J7, J8, J9, J10,
     *                  J11
C     .. Save statement ..
      SAVE              /AD02NM/, /BD02NM/, /CD03PC/, /FD03PC/,
     *                  /GD03PC/, /HD03PC/, /MD03PC/, /PD03PC/,
     *                  /QD03PC/, /RD03PC/, /SD03PC/
C     .. Executable Statements ..
C
C     Check for step failure ..
C
      IF (IMON.LT.0) RETURN
      IF ( .NOT. REMESH) GO TO 20
      IF (NRMESH.NE.0) THEN
         ICOUNT = ICOUNT + 1
         IF (ICOUNT.EQ.NRMESH .OR. (ICOUNT+NRMESH).EQ.0) GO TO 40
      ELSE IF (NRMESH.EQ.0) THEN
C        *(VP)**  IF(ABS((T-TOUT)/TOUT).LT.UROUND) GO TO 17
         IF (T.GE.TOUT .AND. ICOUNT.EQ.0) THEN
            ICOUNT = 1
            GO TO 60
         END IF
      END IF
C     Normal exit to continue integration.
   20 IMON = 1
      RETURN
C   --------------------------------------------------------------------
C     Remeshing followed by interpolation
C     Workspace organisation - The parts of the workspace concerned with
C     **********************   remeshing are organised as follows:
C     RESWK(I8)  ;  The start of the current mesh.
C     RESWK(I9)  ;  The start of the PDE fluxes.
C    RESWK(I10) ;  The start of the new mesh found by D03PRJ. Array XNP.
C     RESWK(I11) ;  The start of the mesh spacings. Array DXP.
C     RESWK(I12) ;  The start of the monitor function array FMON.
C     RESWK(I13) ;  The start of the array CHISTR.
C                 From RESWK(I11) onwards the workspace is also used in
C                 the call to the cubic spline routine D03PRT.
C     RESWK(I14) ;
C     TO          Used for remeshing with material interfaces.
C     RESWK(I16)
C  ---------------------------------------------------------------------
   40 CONTINUE
      ICOUNT = 0
   60 IF (ITRACE.GT.0 .OR. IPMINF.GT.0) THEN
         WRITE (REC,FMT=99999) T
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
      END IF
      NOLD = NPTS*NPDE
      NIP = NPTS
C    If NOP=0 then NOP calculated by D03PRJ else NOP set = NIP by D03PRJ
C     In which case CONST criterion may be over or under satisfied
      NOP = 0
      IF (NOPUSR.NE.0) NOP = NIP
      IFAIL = 0
C
      MIDNPT = NPTS + 1
      CALL MONFKB(T,NIP,NPDE,RESWK(I8),Y,RESWK(I12))
      CALL D03PRH(T,M,NPDE,NIP,RESWK(I8),MAXNPT,CONST,STPRAT,IPMINF,
     *            RESWK(I12),NOP,RESWK(I10),RESWK(I11),RESWK(I13),IFAIL,
     *            RESWK(I14),NXFIX,IXFIX)
      IF (IFAIL.NE.0) THEN
         IMON = -2
         RETURN
      END IF
C
C     If NRMESH < 0 option test the new mesh to see if it is needed.
C
      IF (RMTYPE.EQ.'REMSET') THEN
         CALL D03PRK(RESWK(I8),NIP,RESWK(I10),NOP,I,DXMESH)
         IF (I.EQ.0) THEN
            IMON = 1
            RETURN
         END IF
      END IF
C     New no. of pts and new value of NEQN
      NPTS = NOP
      N = NPTS*NPDE
      NVST = NPTS*NPDE + 1
      NEQN = N + NV
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99998) NOP
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(RESWK(I10),NOP,29)
      END IF
C ----------------------------------------------------------------------
C  Interpolate Y onto new mesh and put new Y from R into Y using D03PRQ.
C-----------------------------------------------------------------------
      CALL D03PRP(RESWK(I10),R,NOP,RESWK(I8),M,Y,NIP,NPDE,IFAIL,
     *            RESWK(I11),IXFIX,NXFIX)
      CALL D03PRQ(Y,R,N,NOLD,NV)
C-----------------------------------------------------------------------
C     Interpolate YDOT (YDOT is overwritten with new YDOT from R.
C-----------------------------------------------------------------------
      CALL D03PRP(RESWK(I10),R,NOP,RESWK(I8),M,YDOT,NIP,NPDE,IFAIL,
     *            RESWK(I11),IXFIX,NXFIX)
      CALL D03PRQ(YDOT,R,N,NOLD,NV)
C     DO 80 K = 1, NINTER
      DO 80 K = 1, 4
C        Allows for THETA method use (as well as Adams/Gear)
C        May need mods for other integrators.
C        Interpolate YSAVE (and overwrite with new YSAVE)
         CALL D03PRP(RESWK(I10),R,NOP,RESWK(I8),M,YSAVE(1,K),NIP,NPDE,
     *               IFAIL,RESWK(I11),IXFIX,NXFIX)
         CALL D03PRQ(YSAVE(1,K),R,N,NOLD,NV)
   80 CONTINUE
      CALL D03PRM(RESWK(I10),NOP,RESWK(I14),NXFIX,IXFIX,NF)
C      UPDATE RESWK WITH XNP(NOP)
      DO 100 I = 1, NOP
         RESWK(I8+I-1) = RESWK(I10+I-1)
  100 CONTINUE
C
      IMON = 3
C     IMON=3 tells integrator to restart with given YSAVE values
      IF (NF.LT.0) IMON = -2
      RETURN
C
99999 FORMAT (' REMESHING PROCESS COMMENCED AT TIME T = ',D12.4)
99998 FORMAT (' ',I3,' NEW MESH POINTS ')
      END
