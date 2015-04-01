      SUBROUTINE E04UDS(CENTRL,INFORM,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                  FDINT,FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,
     *                  C1,C2,CJAC,CJACU,GRAD,GRADU,HFORWD,HCNTRL,X,
     *                  IUSER,USER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1095 (JUL 1993).
C     MARK 17 REVISED. IER-1615 (JUN 1995).
C
C     ******************************************************************
C     E04UDS evaluates any missing gradients.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written 3-July-1986.
C     This version of E04UDS dated 14-Sep-92.
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
      DOUBLE PRECISION  BIGBND, CDINT, FDINT, FDNORM, OBJF
      INTEGER           INFORM, LDCJ, LDCJU, N, NCNLN
      LOGICAL           CENTRL
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), C2(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), GRAD(N), GRADU(N), HCNTRL(N),
     *                  HFORWD(N), USER(*), X(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      INTEGER           LFDSET, LVLDIF, NCDIFF, NFDIFF
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, DELTA, OBJF1, OBJF2, STEPBL,
     *                  STEPBU, XJ
      INTEGER           I, J, MODE, NCOLJ, NSTATE
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
      DO 80 J = 1, N
         XJ = X(J)
         NCOLJ = 0
         IF (NCDIFF.GT.0) THEN
            DO 20 I = 1, NCNLN
               IF (CJACU(I,J).EQ.RDUMMY) THEN
                  NEEDC(I) = 1
                  NCOLJ = NCOLJ + 1
               ELSE
                  NEEDC(I) = 0
               END IF
   20       CONTINUE
         END IF
C
         IF (NCOLJ.GT.0 .OR. GRADU(J).EQ.RDUMMY) THEN
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
            IF (NCOLJ.GT.0) THEN
               CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,NSTATE,
     *                     IUSER,USER)
               IF (MODE.LT.0) GO TO 100
            END IF
C
            IF (GRADU(J).EQ.RDUMMY) THEN
               CALL OBJFUN(MODE,N,X,OBJF1,GRADU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 100
            END IF
C
            IF (CENTRL) THEN
C              ---------------------------------------------------------
C              Central differences.
C              ---------------------------------------------------------
               X(J) = XJ + DELTA + DELTA
C
               IF (NCOLJ.GT.0) THEN
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C2,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 100
C
                  DO 40 I = 1, NCNLN
                     IF (NEEDC(I).EQ.1) CJAC(I,J) = (FOUR*C1(I)
     *                   -THREE*C(I)-C2(I))/(DELTA+DELTA)
   40             CONTINUE
               END IF
C
               IF (GRADU(J).EQ.RDUMMY) THEN
                  CALL OBJFUN(MODE,N,X,OBJF2,GRADU,NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 100
C
                  GRAD(J) = (FOUR*OBJF1-THREE*OBJF-OBJF2)/(DELTA+DELTA)
C
               END IF
            ELSE
C              ---------------------------------------------------------
C              Forward Differences.
C              ---------------------------------------------------------
               IF (NCOLJ.GT.0) THEN
                  DO 60 I = 1, NCNLN
                     IF (NEEDC(I).EQ.1) CJAC(I,J) = (C1(I)-C(I))/DELTA
   60             CONTINUE
               END IF
C
               IF (GRADU(J).EQ.RDUMMY) GRAD(J) = (OBJF1-OBJF)/DELTA
C
            END IF
         END IF
         X(J) = XJ
C
   80 CONTINUE
C
      RETURN
C
  100 INFORM = MODE
      RETURN
C
C     End of  E04UDS. (NPFD)
C
      END
