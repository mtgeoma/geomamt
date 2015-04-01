      SUBROUTINE E04UCV(MODE,LSUMRY,NEEDFD,UNITQ,LVLDER,LDCJ,LDCJU,LDZY,
     *                  N,NCNLN,NFREE,NZ,BIGBND,CVNORM,EPSRF,FDNORM,
     *                  OBJF,KX,NEEDC,CONFUN,OBJFUN,BL,BU,C,C1,C2,CJAC,
     *                  CJACU,GRAD,GRADU,GQ,HFORWD,HCNTRL,ZY,X,WORK,
     *                  IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C***********************************************************************
C     E04UCV   evaluates any missing gradients and the transformed
C     gradient of the objective function.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written 28-September-1985.
C     This version of  E04UCV    dated  18-January-1986.
C***********************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D0)
      DOUBLE PRECISION  ZERO, HALF
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0)
      DOUBLE PRECISION  THREE, FOUR
      PARAMETER         (THREE=3.0D+0,FOUR=4.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CVNORM, EPSRF, FDNORM, OBJF
      INTEGER           LDCJ, LDCJU, LDZY, LVLDER, MODE, N, NCNLN,
     *                  NFREE, NZ
      LOGICAL           NEEDFD, UNITQ
      CHARACTER*3       LSUMRY
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), C2(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), GQ(N), GRAD(N), GRADU(N),
     *                  HCNTRL(N), HFORWD(N), USER(*), WORK(N), X(N),
     *                  ZY(LDZY,LDZY)
      INTEGER           IUSER(*), KX(N), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      INTEGER           LFDSET, LVLDIF, NCDIFF, NFDIFF, NOUT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CNORM, DELTA, GZNORM, OBJF1,
     *                  OBJF2, STEPBL, STEPBU, XJ
      INTEGER           I, J, NFOUND, NSTATE
      LOGICAL           CENTRL, GOODGQ
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          E04NBW, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
      NSTATE = 0
      MODE = 0
C
C     +    REPEAT
   20 CENTRL = LVLDIF .EQ. 2
C
      IF (NEEDFD) THEN
C           ============================================================
C           Compute any missing constraint gradients.
C           ============================================================
         BIGLOW = -BIGBND
         BIGUPP = BIGBND
C
         DO 80 J = 1, N
            XJ = X(J)
C
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
            IF (CENTRL) THEN
               DELTA = HCNTRL(J)
            ELSE
               DELTA = HFORWD(J)
            END IF
C
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) DELTA = -DELTA
C
C              Only compute the constraint values for the indices
C              corresponding to missing Jacobian elements.
C
            NFOUND = 0
            IF (LVLDER.EQ.0 .OR. LVLDER.EQ.1) THEN
               DO 40 I = 1, NCNLN
                  IF (CJACU(I,J).EQ.RDUMMY) THEN
                     NEEDC(I) = 1
                     NFOUND = NFOUND + 1
                  ELSE
                     NEEDC(I) = 0
                  END IF
   40          CONTINUE
            END IF
C
            IF (NFOUND.GT.0) THEN
               X(J) = XJ + DELTA
               CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,NSTATE,
     *                     IUSER,USER)
               IF (MODE.LT.0) GO TO 100
C
               IF (CENTRL) THEN
                  X(J) = XJ + DELTA + DELTA
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C2,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 100
               END IF
C
               DO 60 I = 1, NCNLN
                  IF (NEEDC(I).EQ.1) THEN
                     IF (CENTRL) THEN
                        CJAC(I,J) = (FOUR*C1(I)-THREE*C(I)-C2(I))
     *                              /(DELTA+DELTA)
                     ELSE
                        CJAC(I,J) = (C1(I)-C(I))/DELTA
                     END IF
                  END IF
   60          CONTINUE
            END IF
C
C              ---------------------------------------------------------
C              Repeat for the objective gradients.
C              ---------------------------------------------------------
            IF (GRADU(J).EQ.RDUMMY) THEN
               X(J) = XJ + DELTA
               CALL OBJFUN(MODE,N,X,OBJF1,GRADU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 100
C
               IF (CENTRL) THEN
                  X(J) = XJ + DELTA + DELTA
                  CALL OBJFUN(MODE,N,X,OBJF2,GRADU,NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 100
C
                  GRAD(J) = (FOUR*OBJF1-THREE*OBJF-OBJF2)/(DELTA+DELTA)
               ELSE
                  GRAD(J) = (OBJF1-OBJF)/DELTA
               END IF
            END IF
            X(J) = XJ
   80    CONTINUE
      END IF
C
C        ---------------------------------------------------------------
C        Install the transformed gradient of the objective.
C        ---------------------------------------------------------------
      CALL DCOPY(N,GRAD,1,GQ,1)
      CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,GQ,ZY,WORK)
C
      IF (NEEDFD) THEN
C
C           If X is close to a K-T point,  switch to central differences
C           and recompute the derivatives.
C
         GZNORM = ZERO
         CNORM = ZERO
         IF (NZ.GT.0) GZNORM = DNRM2(NZ,GQ,1)
         IF (NCNLN.GT.0) CNORM = DNRM2(NCNLN,C,1)
C
         GOODGQ = GZNORM .GT. ABS(OBJF)*EPSRF/FDNORM .OR. CVNORM .GT.
     *            CNORM*EPSRF/FDNORM
C
         IF ( .NOT. GOODGQ .AND. LVLDIF.EQ.1) LVLDIF = 2
C
      ELSE
         GOODGQ = .TRUE.
      END IF
C
C     +    UNTIL     (GOODGQ  .OR.  CENTRL)
      IF ( .NOT. (GOODGQ .OR. CENTRL)) GO TO 20
C
      IF (CENTRL) LSUMRY(3:3) = 'Central differences'
C
  100 RETURN
C
C     End of  E04UCV. (NPGQ)
C
      END
