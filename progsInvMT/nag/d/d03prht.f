      SUBROUTINE D03PRH(TIME,M,NPDE,NIP,XOP,MAXNPT,CONST,STPRAT,IPMINF,
     *                  FMON,NOP,XNP,DXW,CHISTR,IFDUM,XFIX,NFIX,IXFIX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     MARK 17 REVISED. IER-1553 (JUN 1995).
C ----------------------------------------------------------------------
C     PADRIV routine from SPRINT.
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONST, STPRAT, TIME
      INTEGER           IFDUM, IPMINF, M, MAXNPT, NFIX, NIP, NOP, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  CHISTR(MAXNPT), DXW(MAXNPT), FMON(NIP), XFIX(*),
     *                  XNP(MAXNPT), XOP(NIP)
      INTEGER           IXFIX(*)
C     .. Scalars in Common ..
      INTEGER           IDEV, ITRACE
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST1
      INTEGER           I, IX, IXN, MANPTI, NIP1, NOP1, NXFXP1
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          D02NNN, D03PRJ, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
C     .. Save statement ..
      SAVE              /AD02NM/
C     .. Executable Statements ..
      IF (NFIX.EQ.0) THEN
         CALL D03PRJ(TIME,M,NPDE,NIP,XOP,MAXNPT,CONST,STPRAT,IPMINF,
     *               FMON,NOP,XNP,DXW,CHISTR,IFDUM)
C VP NOV94
         IF (IFDUM.NE.0) RETURN
      ELSE
         IXN = 1
         NXFXP1 = NFIX + 1
C        WRITE(IDEV,13)(IXFIX(I),I=1,NFIX)
C        13      FORMAT(' FIXED POINT INDICES ARE ',I3)
         DO 20 I = 1, NXFXP1
            IF (I.EQ.1) THEN
               NIP1 = IXFIX(1) + 0.1D0
               IX = 1
            ELSE IF (I.EQ.NXFXP1) THEN
               NIP1 = NIP - IXFIX(NFIX) + 1.1D0
               IX = IXFIX(NFIX)
            ELSE
               NIP1 = IXFIX(I) - IXFIX(I-1) + 1.1D0
               IX = IXFIX(I-1)
            END IF
C
C           Set MANPTI so that new mesh points are in the same ratio
C           as old ones
C
            MANPTI = MAX(MAXNPT*NIP1/NIP,NIP1)
            CONST1 = CONST*(NIP-1.D0)/(NIP1-1.D0)
            NOP1 = 1
            IF (NOP.EQ.0) NOP1 = 0
            IF (IPMINF.GE.1) THEN
               CALL X04ABF(0,IDEV)
               WRITE (REC,FMT=99999)
               CALL X04BAF(IDEV,REC)
               CALL D02NNN(XOP(IX),NIP1,29)
               WRITE (REC,FMT=99998) NIP1
               CALL X04BAF(IDEV,REC)
               WRITE (REC,FMT=99997) IX
               CALL X04BAF(IDEV,REC)
            END IF
            CALL D03PRJ(TIME,M,NPDE,NIP1,XOP(IX),MANPTI,CONST1,STPRAT,
     *                  IPMINF,FMON(IX),NOP1,XNP(IXN),DXW(IXN),
     *                  CHISTR(IXN),IFDUM)
C VP NOV94
            IF (IFDUM.NE.0) RETURN
C
C           Ensure interface points are O.K.
C
            XNP(IXN) = XOP(IX)
            XNP(IXN+NOP1-1) = XOP(IX+NIP1-1)
            IF (IPMINF.GE.1) THEN
               WRITE (REC,FMT=99996) NOP1
               CALL X04BAF(IDEV,REC)
               WRITE (REC,FMT=99995)
               CALL X04BAF(IDEV,REC)
               CALL D02NNN(XNP(IXN),NOP1,29)
            END IF
            IXN = IXN + NOP1 - 1
   20    CONTINUE
         NOP = IXN
         IF (IPMINF.GE.1) THEN
            WRITE (REC,FMT=99994) NOP
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
      RETURN
C
99999 FORMAT (' OLD MESH')
99998 FORMAT (' D03PRJ CALLED WITH',I4,' POINTS')
99997 FORMAT (' START POSITION OF MESH SEGMENT = ',I4)
99996 FORMAT (' NO OF OUTPUT PTS FROM D03PRJ = ',I5)
99995 FORMAT (' NEW MESH')
99994 FORMAT (' FINAL NO OF OUTPUT POINTS FROM D03PRH =',I5)
      END
