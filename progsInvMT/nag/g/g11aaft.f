      SUBROUTINE G11AAF(NROW,NCOL,NOBST,LDT,EXPT,CHIST,PROB,CHI,G,DF,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G11AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHI, DF, G, PROB
      INTEGER           IFAIL, LDT, NCOL, NROW
C     .. Array Arguments ..
      DOUBLE PRECISION  CHIST(LDT,NCOL), EXPT(LDT,NCOL)
      INTEGER           NOBST(LDT,NCOL)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, C, CC, E, O, PE, PG, PLE
      INTEGER           I, IERROR, IFAULT, J, L, M, MC1, MC2, MR1, MR2,
     *                  N, NREC
      LOGICAL           EXACT, WARN
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01BLF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, LOG
C     .. Executable Statements ..
C
C     initialize and check inputs
C
      NREC = 1
      WARN = .FALSE.
      EXACT = .FALSE.
      IERROR = 1
      IF (NROW.LT.2) THEN
         WRITE (P01REC,FMT=99999) NROW
      ELSE IF (NCOL.LT.2) THEN
         WRITE (P01REC,FMT=99998) NCOL
      ELSE IF (LDT.LT.NROW) THEN
         WRITE (P01REC,FMT=99997) LDT, NROW
      ELSE
         IERROR = 0
C
C        Compute  row margin and total and check observations
C
         N = 0
         DO 40 I = 1, NROW
            M = 0
            DO 20 J = 1, NCOL
               L = NOBST(I,J)
               IF (L.LT.0) THEN
                  IERROR = 2
                  WRITE (P01REC,FMT=99996)
                  GO TO 120
               END IF
               M = M + L
   20       CONTINUE
            EXPT(I,NCOL) = M
            N = N + M
   40    CONTINUE
         IF (N.EQ.0) THEN
            IERROR = 2
            WRITE (P01REC,FMT=99995)
            GO TO 120
         END IF
C
C        Compute expected frequencies, chi-square contributions
C        chi-square statistic and likelihood ratio
C
         CHI = 0.0D0
         G = 0.0D0
         DO 100 I = 1, NCOL
C
C           Compute column margin
C
            M = 0
            DO 60 J = 1, NROW
               M = M + NOBST(J,I)
   60       CONTINUE
            C = DBLE(M)/DBLE(N)
            DO 80 J = 1, NROW
               E = C*EXPT(J,NCOL)
               EXPT(J,I) = E
               IF (E.LE.0.0D0) THEN
                  CHIST(J,I) = 0.0D0
                  WARN = .TRUE.
               ELSE
                  O = DBLE(NOBST(J,I))
                  CC = (O-E)*(O-E)/E
                  CHIST(J,I) = CC
                  CHI = CHI + CC
                  IF (O.GT.0.0D0) G = G + O*LOG(O/E)
                  IF (E.LE.0.5D0) WARN = .TRUE.
               END IF
   80       CONTINUE
  100    CONTINUE
         G = 2.0D0*G
C
C        Compute probability
C
         IF (NROW.EQ.2 .AND. NCOL.EQ.2) THEN
C
C           For 2 X 2 use Yates correction
C
            A = NOBST(1,1)*NOBST(2,2) - NOBST(1,2)*NOBST(2,1)
            A = ABS(A) - 0.5D0*DBLE(N)
            MC1 = NOBST(1,1) + NOBST(2,1)
            MC2 = NOBST(1,2) + NOBST(2,2)
            MR1 = NOBST(1,1) + NOBST(1,2)
            MR2 = NOBST(2,1) + NOBST(2,2)
            IF (MR1.EQ.0 .OR. MR2.EQ.0 .OR. MC1.EQ.0 .OR. MC2.EQ.0) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99993)
               GO TO 120
            END IF
            CHI = (DBLE(N)/DBLE(MR1*MR2))*((A*A)/DBLE(MC1*MC2))
            DF = 1.0D0
            IF (N.GT.40) THEN
               IFAULT = 1
               PROB = G01ECF('U',CHI,DF,IFAULT)
            ELSE
C
C                 For 2 X 2 with N small use Fisher exact
C
               IFAULT = 1
               CALL G01BLF(N,MC1,MR1,NOBST(1,1),PLE,PG,PE,IFAULT)
               PROB = PG + PE
               IF (PROB.GT.PLE) PROB = PLE
               IF (PROB.LT.0.5D0) THEN
                  PROB = 2.0D0*PROB
               ELSE
                  PROB = 1.0D0
               END IF
               EXACT = .TRUE.
            END IF
         ELSE
C
C           Compute chi-square probability for general case
C
            IFAULT = 1
            DF = (NCOL-1)*(NROW-1)
            PROB = G01ECF('U',CHI,DF,IFAULT)
         END IF
         IF (WARN .AND. ( .NOT. EXACT)) THEN
C
C           Set warning for low cell frequency
C
            IERROR = 4
            WRITE (P01REC,FMT=99994)
         END IF
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
99999 FORMAT (' ** On entry, NROW.lt.2: NROW = ',I16)
99998 FORMAT (' ** On entry, NCOL.lt.2: NCOL = ',I16)
99997 FORMAT (' ** On entry, LDT.lt.NROW: LDT = ',I16,' NROW = ',I16)
99996 FORMAT (' ** On entry, at least one element of NOBST.lt.0')
99995 FORMAT (' ** On entry, all elements of NOBST = 0')
99994 FORMAT (' ** Warning: a cell has an expected frequency .le. 0.5')
99993 FORMAT (' ** On entry, a 2 x 2 table is degenerate')
      END
