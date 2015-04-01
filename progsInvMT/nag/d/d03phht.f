      SUBROUTINE D03PHH(PHFPDE,PDEPCF,T,M,U,NPDE,NPTS,RESWK,NRESWK,V,NV,
     *                  MAXNPT,IRES,XI,NXI,XFIX,IXFIX,NXFIX,NF,UVINIT,
     *                  SPDEF1,MONFFD,MONFKB)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1993.
C     MARK 17 REVISED. IER-1549 (JUN 1995).
C
C  ---------------------------------------------------------------------
C     Routine to provide initial mesh by three calls to D03PRH.
C     (Finite Differences)
C     (INMESH routine from SPRINT)
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IRES, M, MAXNPT, NF, NPDE, NPTS, NRESWK, NV,
     *                  NXFIX, NXI
C     .. Array Arguments ..
      DOUBLE PRECISION  RESWK(NRESWK), U(NPDE,*), V(*), XFIX(*), XI(*)
      INTEGER           IXFIX(*)
C     .. Subroutine Arguments ..
      EXTERNAL          MONFFD, MONFKB, PDEPCF, PHFPDE, SPDEF1, UVINIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  CONST, DYMESH, STPRAT, TOUT
      INTEGER           I1, I10, I11, I12, I13, I14, I15, I16, I2, I3,
     *                  I4, I5, I6, I7, I8, I9, ICOUNT, IDEV, IIFLAG,
     *                  IPMINF, ITRACE, NOPUSR, NRMESH
      LOGICAL           REMESH
      CHARACTER*6       PDCODE, RMTYPE
C     .. Local Scalars ..
      DOUBLE PRECISION  XP
      INTEGER           I, IFAIL, IONE, ISTATE, ITWO, J, K, K1, K2, K3,
     *                  K4, MIDNPT, NIP, NIPP1, NOP
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  XXP(1)
C     .. External Subroutines ..
      EXTERNAL          D02NNN, D03PRH, D03PRK, D03PRM, D03PZW, X04ABF,
     *                  X04BAF
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /GD03PC/CONST, STPRAT, IPMINF, NOPUSR
      COMMON            /HD03PC/TOUT, REMESH, NRMESH, ICOUNT
      COMMON            /MD03PC/DYMESH
      COMMON            /QD03PC/I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     *                  I11, I12, I13, I14, I15, I16
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /AD02NM/, /FD03PC/, /GD03PC/, /HD03PC/,
     *                  /MD03PC/, /QD03PC/, /XD03PC/
C     .. Executable Statements ..
C
C     If NOP=0 then NOP calculated by D03PRJ (called by D03PRH) else
C     NOP set equal to NIP by D03PRJ and so const criterion may be
C     over or under satisfied.
C
      CALL D03PRM(RESWK(I8),NPTS,XFIX,NXFIX,IXFIX,NF)
      IF (NF.LT.0) RETURN
      K1 = NPDE + 1
      K2 = K1 + NPDE
      K3 = K2 + NPDE*NPDE
      K4 = K3 + NPDE
      NIP = NPTS
      DO 100 K = 1, 3
C        (VP)  NOP=0
         NOP = NIP
C VP NOV94
         IFAIL = 0
         ISTATE = 0
         DO 20 I = 1, NIP
            RESWK(I10+I-1) = RESWK(I8+I-1)
   20    CONTINUE
         IF (K.EQ.1 .AND. ITRACE.GE.1) THEN
            CALL X04ABF(0,IDEV)
            WRITE (REC,FMT=99999)
            CALL X04BAF(IDEV,REC)
            CALL D02NNN(RESWK(I10),NIP,29)
         END IF
C ----------------------------------------------------------------------
C    Formation of first derivatives and then the flux by calling D03PZW.
C ----------------------------------------------------------------------
         IONE = 1
         ITWO = 2
         NIPP1 = NIP + 1
         DO 60 I = 1, NIPP1
            IF (I.EQ.1) THEN
               XXP(1) = RESWK(I10)
               XP = XXP(1)
            ELSE IF (I.EQ.NIPP1) THEN
               XXP(1) = RESWK(I10+NIP-1)
               XP = XXP(1)
            ELSE
               XXP(1) = (RESWK(I10+I)+RESWK(I10+I-1))*0.5D0
               XP = XXP(1)
            END IF
            CALL D03PZW(XXP,RESWK,IONE,RESWK(I10),M,U,NIP,NPDE,ITWO,
     *                  IFAIL)
            IF (PDCODE.EQ.'SPKEEL') THEN
               CALL SPDEF1(PHFPDE,PDEPCF,T,XP,NPDE,RESWK,RESWK(K1),
     *                     RESWK(K2),RESWK(K3),RESWK(K4),NV,V,RESWK(I7),
     *                     IRES,IIFLAG)
               IF (IRES.NE.1) RETURN
               DO 40 J = 1, NPDE
                  RESWK(I9+(I-1)*NPDE+J-1) = RESWK(K4+J-1)
   40          CONTINUE
            END IF
   60    CONTINUE
         MIDNPT = NIP + 1
         IF (PDCODE.EQ.'SPKEEL') THEN
            CALL MONFFD(T,NIP,NPDE,RESWK(I10),U,RESWK(I9),RESWK(I12))
         ELSE
            CALL MONFKB(T,NIP,NPDE,RESWK(I10),U,RESWK(I12))
         END IF
         CALL D03PRH(T,M,NPDE,NIP,RESWK(I8),MAXNPT,CONST,STPRAT,IPMINF,
     *               RESWK(I12),NOP,RESWK(I10),RESWK(I11),RESWK(I13),
     *               IFAIL,XFIX,NXFIX,IXFIX)
C VP NOV94
         IF (IFAIL.NE.0) THEN
            NF = -2
            RETURN
         END IF
C
         IF (RMTYPE.EQ.'REMSET') THEN
C           Test if new mesh really needed
            CALL D03PRK(RESWK(I8),NIP,RESWK(I10),NOP,I,DYMESH)
            IF (I.EQ.0) THEN
C               New mesh not needed
               RETURN
            ELSE IF (ITRACE.GT.0 .OR. IPMINF.GT.0) THEN
C               New mesh is needed
               IF (ITRACE.GE.1) THEN
                  WRITE (REC,FMT=99998)
                  CALL X04BAF(IDEV,REC)
               END IF
            END IF
         ELSE IF (NRMESH.EQ.0 .AND. T.LT.TOUT) THEN
C             New mesh not needed
            RETURN
         END IF
         IF (ITRACE.GE.1) THEN
            WRITE (REC,FMT=99997)
            CALL X04BAF(IDEV,REC)
            CALL D02NNN(RESWK(I10),NOP,29)
         END IF
C  ---------------------------------------------------------------------
C        Find the correct initial values at new mesh points.
C  ---------------------------------------------------------------------
         DO 80 I = 1, NOP
            RESWK(I8+I-1) = RESWK(I10+I-1)
   80    CONTINUE
C
         CALL UVINIT(NPDE,NOP,NXI,RESWK(I10),XI,U,NV,V)
         CALL D03PRM(RESWK(I8),NOP,XFIX,NXFIX,IXFIX,NF)
         IF (NF.LT.0) RETURN
         NIP = NOP
         NPTS = NOP
  100 CONTINUE
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99996)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(RESWK(I10),NIP,29)
      END IF
      RETURN
C
99999 FORMAT ('  INITIAL MESH POINTS ')
99998 FORMAT ('  NEW MESH ADOPTED  ')
99997 FORMAT ('  NEW MESH POINTS  ')
99996 FORMAT ('  FINAL MESH POINTS ')
      END
