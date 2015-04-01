      SUBROUTINE E04VDQ(NFREE,NROWA,NROWJ,N,NCLIN,NCNLN,NCTOTL,BIGBND,
     *                  NAMED,NAMES,LENNAM,NACTIV,ISTATE,KACTIV,A,BL,BU,
     *                  C,CLAMDA,RLAMDA,X)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. IER-620 (APR 1988).
C
C *********************************************************************
C     E04VDQ  EXPANDS THE LAGRANGE MULTIPLIERS INTO  CLAMDA.
C     IF  MSG .GE. 10  OR  MSG .EQ. 1,  E04VDQ  THEN PRINTS  X, A*X,
C     C(X), THEIR BOUNDS,  THE MULTIPLIERS, AND THE RESIDUALS
C     (DISTANCE TO THE NEAREST BOUND).
C     E04VDQ  IS CALLED BY  E04MBY, E04NAX, LCCORE AND E04VCZ  JUST
C     BEFORE THEY EXIT.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF MARCH 1982. REV. OCT. 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           LENNAM, N, NACTIV, NCLIN, NCNLN, NCTOTL, NFREE,
     *                  NROWA, NROWJ
      LOGICAL           NAMED
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), BL(NCTOTL), BU(NCTOTL), C(NROWJ),
     *                  CLAMDA(NCTOTL), RLAMDA(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), NAMES(4,LENNAM)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RES, RES2, V, WLAM, ZERO
      INTEGER           IP, IS, J, JJ, K, LROWA, NFIXED, NLAM, NPLIN
      CHARACTER*1       ID3
      CHARACTER*2       LS
C     .. Local Arrays ..
      CHARACTER*1       ID(3)
      CHARACTER*2       LSTATE(7)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FBF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
C     .. Data statements ..
      DATA              ID(1), ID(2), ID(3)/'V', 'L', 'N'/
      DATA              LSTATE(1), LSTATE(2)/'--', '++'/
      DATA              LSTATE(3), LSTATE(4)/'FR', 'LL'/
      DATA              LSTATE(5), LSTATE(6)/'UL', 'EQ'/
      DATA              LSTATE(7)/'TB'/
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
C
      NPLIN = N + NCLIN
      LROWA = NROWA*(N-1) + 1
C
C     EXPAND BOUND, LINEAR AND NONLINEAR MULTIPLIERS INTO  CLAMDA.
C
      CALL F06FBF(NCTOTL,ZERO,CLAMDA,1)
      NFIXED = N - NFREE
      NLAM = NACTIV + NFIXED
      IF (NLAM.EQ.0) GO TO 40
C
      DO 20 K = 1, NLAM
         J = KACTIV(K)
         IF (K.LE.NACTIV) J = J + N
         CLAMDA(J) = RLAMDA(K)
   20 CONTINUE
C
   40 IF (MSG.LT.10 .AND. MSG.NE.1) RETURN
C
      WRITE (REC,FMT=99999)
      DO 60 J = 1, 4
         CALL X04BAF(NOUT,REC(J))
   60 CONTINUE
      ID3 = ID(1)
C
      DO 260 J = 1, NCTOTL
         B1 = BL(J)
         B2 = BU(J)
         WLAM = CLAMDA(J)
         IS = ISTATE(J)
         LS = LSTATE(IS+3)
         IF (J.LE.N) GO TO 80
         IF (J.LE.NPLIN) GO TO 100
         GO TO 160
C
C
C        SECTION 1 -- THE VARIABLES  X.
C        ------------------------------
   80    K = J
         V = X(J)
         GO TO 220
C
C
C        SECTION 2 -- THE LINEAR CONSTRAINTS  A*X.
C        -----------------------------------------
  100    IF (J.NE.N+1) GO TO 140
         WRITE (REC,FMT=99998)
         DO 120 JJ = 1, 4
            CALL X04BAF(NOUT,REC(JJ))
  120    CONTINUE
         ID3 = ID(2)
C
  140    K = J - N
         V = DDOT(N,A(K,1),NROWA,X,1)
         GO TO 220
C
C
C        SECTION 3 -- THE NONLINEAR CONSTRAINTS  C(X).
C        ---------------------------------------------
C
  160    IF (NCNLN.LE.0) GO TO 260
         IF (J.NE.NPLIN+1) GO TO 200
         WRITE (REC,FMT=99997)
         DO 180 JJ = 1, 4
            CALL X04BAF(NOUT,REC(JJ))
  180    CONTINUE
         ID3 = ID(3)
C
  200    K = J - NPLIN
         V = C(K)
C
C
C        PRINT A LINE FOR THE J-TH VARIABLE OR CONSTRAINT.
C        -------------------------------------------------
  220    RES = V - B1
         RES2 = B2 - V
         IF (ABS(RES).GT.ABS(RES2)) RES = RES2
         IP = 1
         IF (B1.LE.(-BIGBND)) IP = 2
         IF (B2.GE.BIGBND) IP = IP + 2
C        THE FOLLOWING UNCONDITIONAL GO TO REPLACES THE COMMENTED
C        OUT CODE. THIS IS TO ENSURE THAT THE NAME OPTION IS NOT
C        USED WITH REVISED NAG OUTPUT FORMATS
         GO TO 240
C        IF (.NOT. NAMED) GO TO 490
C
C        DO 450 L = 1, 4
C           ID4(L) = NAMES(L,J)
C        450    CONTINUE
C        IF (IP .EQ. 1) WRITE (NOUT, 2100) ID4,    LS, V,B1,B2,WLAM,RES
C        IF (IP .EQ. 2) WRITE (NOUT, 2200) ID4,    LS, V,   B2,WLAM,RES
C        IF (IP .EQ. 3) WRITE (NOUT, 2300) ID4,    LS, V,B1,   WLAM,RES
C        IF (IP .EQ. 4) WRITE (NOUT, 2400) ID4,    LS, V,      WLAM,RES
C        GO TO 500
C
  240    IF (IP.EQ.1) WRITE (REC,FMT=99996) ID3, K, LS, V, B1, B2, WLAM,
     *       RES
         IF (IP.EQ.2) WRITE (REC,FMT=99995) ID3, K, LS, V, B2, WLAM, RES
         IF (IP.EQ.3) WRITE (REC,FMT=99994) ID3, K, LS, V, B1, WLAM, RES
         IF (IP.EQ.4) WRITE (REC,FMT=99993) ID3, K, LS, V, WLAM, RES
         IF (IP.GE.1 .AND. IP.LE.4) CALL X04BAF(NOUT,REC(1))
  260 CONTINUE
C
      RETURN
C
C2100 FORMAT(1X, 4A2, 10X, A2, 3G16.7, G16.7, G16.4)
C2200 FORMAT(1X, 4A2, 10X, A2, G16.7, 5X, 5H NONE, 6X, G16.7,
C     *   G16.7, G16.4)
C2300 FORMAT(1X, 4A2, 10X, A2, 2G16.7, 5X, 5H NONE, 6X, G16.7, G16.4)
C     2400 FORMAT(1X, 4A2, 10X, A2,  G16.7, 5X, 5H NONE, 11X, 5H NONE,
C     *   6X, G16.7, G16.4)
C
C     END OF E04VDQ  ( PRTSOL )
99999 FORMAT (//' VARBL STATE',4X,' VALUE',5X,' LOWER BOUND',3X,' UPPE',
     *  'R BOUND    LAGR MULT   RESIDUAL',/)
99998 FORMAT (//' LNCON STATE',4X,' VALUE',5X,' LOWER BOUND',3X,' UPPE',
     *  'R BOUND    LAGR MULT   RESIDUAL',/)
99997 FORMAT (//' NLCON STATE',4X,' VALUE',5X,' LOWER BOUND',3X,' UPPE',
     *  'R BOUND    LAGR MULT   RESIDUAL',/)
99996 FORMAT (1X,1A1,I3,4X,A2,3G15.7,2G12.4)
99995 FORMAT (1X,1A1,I3,4X,A2,G15.7,5X,'NONE',6X,G15.7,2G12.4)
99994 FORMAT (1X,1A1,I3,4X,A2,2G15.7,5X,'NONE',6X,2G12.4)
99993 FORMAT (1X,1A1,I3,4X,A2,G15.7,5X,'NONE',11X,'NONE',6X,2G12.4)
      END
