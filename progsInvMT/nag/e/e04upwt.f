      SUBROUTINE E04UPW(CENTRL,INFORM,LDCJ,LDCJU,LDFJ,LDFJU,M,N,NCNLN,
     *                  BIGBND,CDINT,FDINT,FDNORM,CONFUN,OBJFUN,NEEDC,
     *                  BL,BU,C,C1,C2,CJAC,CJACU,F,F1,F2,FJAC,FJACU,
     *                  HFORWD,HCNTRL,X,IUSER,USER)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1102 (JUL 1993).
C     MARK 17 REVISED. IER-1622 (JUN 1995).
C
C     ******************************************************************
C     E04UPW evaluates any missing gradients.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version based on NPFD written 3-July-1986.
C     This version of E04UPW dated 11-May-1988.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      DOUBLE PRECISION  THREE, FOUR
      PARAMETER         (THREE=3.0D+0,FOUR=4.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CDINT, FDINT, FDNORM
      INTEGER           INFORM, LDCJ, LDCJU, LDFJ, LDFJU, M, N, NCNLN
      LOGICAL           CENTRL
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), C2(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), F(M), F1(M), F2(M),
     *                  FJAC(LDFJ,*), FJACU(LDFJU,*), HCNTRL(N),
     *                  HFORWD(N), USER(*), X(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      INTEGER           LFDSET, LVLDIF, NCDIFF, NFDIFF
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, DELTA, STEPBL, STEPBU, XJ
      INTEGER           I, J, MODE, NCMISS, NFMISS, NSTATE
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
C     .. Executable Statements ..
C
      INFORM = 0
C
C     ==================================================================
C     Use the pre-assigned difference intervals to approximate the
C     derivatives.
C     ==================================================================
C     Use either the same interval for each element (LFDSET = 1),
C     or the intervals already in HFORWD or HCNTRL (LFDSET = 0 or 2).
C
      NSTATE = 0
      MODE = 0
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
      FDNORM = ZERO
C
      DO 140 J = 1, N
C
         XJ = X(J)
C
         NCMISS = 0
         IF (NCDIFF.GT.0) THEN
            DO 20 I = 1, NCNLN
               IF (CJACU(I,J).EQ.RDUMMY) THEN
                  NEEDC(I) = 1
                  NCMISS = NCMISS + 1
               ELSE
                  NEEDC(I) = 0
               END IF
   20       CONTINUE
         END IF
C
         NFMISS = 0
         IF (NFDIFF.GT.0) THEN
            DO 40 I = 1, M
               IF (FJACU(I,J).EQ.RDUMMY) NFMISS = NFMISS + 1
   40       CONTINUE
         END IF
C
         IF (NCMISS.GT.0 .OR. NFMISS.GT.0) THEN
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
            IF (CENTRL) THEN
               IF (LFDSET.EQ.1) THEN
                  DELTA = CDINT
               ELSE
                  DELTA = HCNTRL(J)
               END IF
            ELSE
               IF (LFDSET.EQ.1) THEN
                  DELTA = FDINT
               ELSE
                  DELTA = HFORWD(J)
               END IF
            END IF
C
            DELTA = DELTA*(ONE+ABS(XJ))
            FDNORM = MAX(FDNORM,DELTA)
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) DELTA = -DELTA
C
            X(J) = XJ + DELTA
C
            IF (NCMISS.GT.0) THEN
               CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,NSTATE,
     *                     IUSER,USER)
               IF (MODE.LT.0) GO TO 160
            END IF
C
            IF (NFMISS.GT.0) THEN
               CALL OBJFUN(MODE,M,N,LDFJU,X,F1,FJACU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 160
            END IF
C
            IF (CENTRL) THEN
C              ---------------------------------------------------------
C              Central differences.
C              ---------------------------------------------------------
               X(J) = XJ + DELTA + DELTA
C
               IF (NCMISS.GT.0) THEN
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C2,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 160
C
                  DO 60 I = 1, NCNLN
                     IF (NEEDC(I).EQ.1) CJAC(I,J) = (FOUR*C1(I)
     *                   -THREE*C(I)-C2(I))/(DELTA+DELTA)
   60             CONTINUE
               END IF
C
               IF (NFMISS.GT.0) THEN
                  CALL OBJFUN(MODE,M,N,LDFJU,X,F2,FJACU,NSTATE,IUSER,
     *                        USER)
                  IF (MODE.LT.0) GO TO 160
C
                  DO 80 I = 1, M
                     IF (FJACU(I,J).EQ.RDUMMY) FJAC(I,J) = (FOUR*F1(I)
     *                   -THREE*F(I)-F2(I))/(DELTA+DELTA)
   80             CONTINUE
               END IF
            ELSE
C              ---------------------------------------------------------
C              Forward Differences.
C              ---------------------------------------------------------
               IF (NCMISS.GT.0) THEN
                  DO 100 I = 1, NCNLN
                     IF (NEEDC(I).EQ.1) CJAC(I,J) = (C1(I)-C(I))/DELTA
  100             CONTINUE
               END IF
C
               IF (NFMISS.GT.0) THEN
                  DO 120 I = 1, M
                     IF (FJACU(I,J).EQ.RDUMMY) FJAC(I,J) = (F1(I)-F(I))
     *                   /DELTA
  120             CONTINUE
               END IF
            END IF
         END IF
         X(J) = XJ
C
  140 CONTINUE
C
      RETURN
C
  160 INFORM = MODE
      RETURN
C
C     End of  E04UPW.  (NLFD)
C
      END
