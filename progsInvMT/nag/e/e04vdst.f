      SUBROUTINE E04VDS(LPROB,N,NCLIN0,NCTOTL,NACTIV,NCOLZ,NFREE,NROWA,
     *                  NROWRT,NCOLRT,JSMLST,KSMLST,SMLLST,ISTATE,
     *                  KACTIV,A,ANORM,QTG,RLAMDA,RT)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11C REVISED. IER-465 (MAR 1985)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C *********************************************************************
C     E04VDS  FIRST COMPUTES THE LAGRANGE MULTIPLIER ESTIMATES FOR THE
C     GIVEN WORKING SET.  IT THEN DETERMINES THE VALUES AND INDICES OF
C     CERTAIN SIGNIFICANT MULTIPLIERS.  IN THIS PROCESS, THE MULTIPLIERS
C     FOR INEQUALITIES AT THEIR UPPER BOUNDS ARE ADJUSTED SO THAT A
C     NEGATIVE MULTIPLIER FOR AN INEQUALITY CONSTRAINT INDICATES
C     NON-OPTIMALITY.  IN THE FOLLOWING, THE TERM MINIMUM REFERS TO THE
C     ORDERING OF NUMBERS ON THE REAL LINE, AND NOT TO THEIR MAGNITUDE.
C
C     SMLLST  IS THE MINIMUM AMONG THE INEQUALITY CONSTRAINTS OF THE
C          (ADJUSTED) MULTIPLIERS SCALED BY THE 2-NORM OF THE
C          ASSOCIATED CONSTRAINT ROW.
C
C     JSMLST  IS THE INDEX OF THE CONSTRAINT CORRESPONDING TO  SMLLST.
C     KSMLST  MARKS ITS POSITION IN  KACTIV.
C
C
C     ON EXIT,  ELEMENTS  1  THRU  NACTIV  OF   RLAMDA  CONTAIN THE
C     (UNADJUSTED) MULTIPLIERS FOR THE GENERAL CONSTRAINTS.  ELEMENTS
C     NACTIV  ONWARDS OF  RLAMDA  CONTAIN THE (UNADJUSTED) MULTIPLIERS
C     FOR THE SIMPLE BOUNDS.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     ORIGINAL VERSION OCTOBER 1982.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SMLLST
      INTEGER           JSMLST, KSMLST, N, NACTIV, NCLIN0, NCOLRT,
     *                  NCOLZ, NCTOTL, NFREE, NROWA, NROWRT
      CHARACTER*2       LPROB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), ANORM(NCLIN0), QTG(N), RLAMDA(N),
     *                  RT(NROWRT,NCOLRT)
      INTEGER           ISTATE(NCTOTL), KACTIV(N)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORMJ, BLAM, FLMAX, ONE, RLAM
      INTEGER           I, IDIAG, IS, J, JGFXD, K, KA, KB, L, L1, L2,
     *                  NFIXED, NLAM
C     .. Local Arrays ..
      CHARACTER*70      REC(3)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F04YAZ, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ONE/1.0D+0/
C     .. Executable Statements ..
C
      FLMAX = WMACH(7)
C
C ---------------------------------------------------------------------
C     FIRST, COMPUTE THE LAGRANGE MULTIPLIERS FOR THE GENERAL
C     CONSTRAINTS IN THE WORKING SET, BY SOLVING
C     T(TRANSPOSE)*RLAMDA = Y(T)*GRAD.
C ---------------------------------------------------------------------
      NFIXED = N - NFREE
      NLAM = NFIXED + NACTIV
      IF (NACTIV.EQ.0) GO TO 20
      CALL DCOPY(NACTIV,QTG(NCOLZ+1),1,RLAMDA,1)
      IDIAG = 1
      CALL F04YAZ(-2,NACTIV,RT(1,NCOLZ+1),NROWRT,RLAMDA,IDIAG)
C
C ---------------------------------------------------------------------
C     NOW SET ELEMENTS NACTIV, NACTIV+1,... OF RLAMDA EQUAL TO THE
C     MULTIPLIERS FOR THE BOUND CONSTRAINTS IN THE WORKING SET.
C ---------------------------------------------------------------------
   20 IF (NFIXED.EQ.0) GO TO 100
      DO 80 L = 1, NFIXED
         KB = NACTIV + L
         J = KACTIV(KB)
         JGFXD = NFREE + L
         BLAM = QTG(JGFXD)
         IF (NACTIV.EQ.0) GO TO 60
         DO 40 KA = 1, NACTIV
            I = KACTIV(KA)
            BLAM = BLAM - A(I,J)*RLAMDA(KA)
   40    CONTINUE
   60    RLAMDA(KB) = BLAM
   80 CONTINUE
C
C ---------------------------------------------------------------------
C     FIND  ALLMAX  AND  SMLLST.
C ---------------------------------------------------------------------
  100 SMLLST = FLMAX
      JSMLST = 0
      KSMLST = 0
      IF (NLAM.EQ.0) GO TO 140
      DO 120 K = 1, NLAM
         J = KACTIV(K)
         IF (K.GT.NACTIV) ANORMJ = ONE
         IF (K.LE.NACTIV) ANORMJ = ANORM(J)
C
         IF (K.LE.NACTIV) J = J + N
         IS = ISTATE(J)
         RLAM = RLAMDA(K)*ANORMJ
C
C        CHANGE THE SIGN OF THE ESTIMATE IF THE CONSTRAINT IS IN THE
C        WORKING SET (OR VIOLATED) AT ITS UPPER BOUND.
C
         IF (IS.EQ.2) RLAM = -RLAM
         IF (IS.EQ.3) RLAM = ABS(RLAM)
         IF (IS.EQ.4) RLAM = -ABS(RLAM)
C
C
C        SKIP IF THIS IS A FIXED VARIABLE OR AN EQUALITY CONSTRAINT.
C
         IF (IS.EQ.3) GO TO 120
C
C        FIND THE SMALLEST MULTIPLIER FOR THE INEQUALITIES.
C
         IF (SMLLST.LE.RLAM) GO TO 120
         SMLLST = RLAM
         JSMLST = J
         KSMLST = K
  120 CONTINUE
C
C ---------------------------------------------------------------------
C     IF REQUIRED, PRINT THE MULTIPLIERS.
C ---------------------------------------------------------------------
  140 IF (MSG.LT.20) GO TO 200
      IF (NACTIV.GT.0) THEN
         WRITE (REC,FMT=99999) LPROB
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 160 I = 1, NACTIV, 4
            WRITE (REC,FMT=99998) (KACTIV(K),RLAMDA(K)
     *        ,K=I,MIN(NACTIV,I+3))
            CALL X04BAF(NOUT,REC(1))
  160    CONTINUE
      END IF
      L1 = NACTIV + 1
      L2 = NLAM
      IF (L1.LE.L2) THEN
         WRITE (REC,FMT=99997) LPROB
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 180 I = L1, L2, 4
            WRITE (REC,FMT=99998) (KACTIV(K),RLAMDA(K),K=I,MIN(L2,I+4))
            CALL X04BAF(NOUT,REC(1))
  180    CONTINUE
      END IF
      IF (MSG.GE.80) THEN
         WRITE (REC,FMT=99996) JSMLST, SMLLST, KSMLST
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
C
  200 RETURN
C
C
C     END OF E04VDS  ( GETLAM )
99999 FORMAT (/' MULTIPLIERS FOR THE ',A2,' CONSTRAINTS...')
99998 FORMAT (4(I5,1P,D11.2))
99997 FORMAT (/' MULTIPLIERS FOR THE ',A2,' BOUND CONSTRAINTS...')
99996 FORMAT (/' //E04VDS//  JSMLST     SMLLST     KSMLST',/' //E04VDS',
     *  '//  ',I6,1P,D11.2,5X,I6)
      END
