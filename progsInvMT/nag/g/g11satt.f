      SUBROUTINE G11SAT(S,N,LP,RL,X,N2,AX,XN,Q,ALPHA,GAMMA,A,C,GPROB,RR,
     *                  CHISQR,ISHOW,IERROR,OB,Y,XL,IA,PHI,P,NROWXR,OBS,
     *                  EXPP,CHI,IDF,VAR,G,IAA,SIGLEV)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-1120 (JUL 1993).
C
C     CALCULATE FACTOR SCORES AND OTHER DIAGNOSTIC INFORMATION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHI, SIGLEV
      INTEGER           IA, IAA, IDF, IERROR, ISHOW, LP, N, N2, NROWXR,
     *                  Q, S
      LOGICAL           CHISQR, GPROB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N2), ALPHA(N2), AX(20), C(N2), EXPP(IA,N2),
     *                  G(2*N2), GAMMA(N2), OBS(IA,N2), P(S),
     *                  PHI(N2,20), RR(N2,20), VAR(IAA,2*N2), XL(S),
     *                  XN(20), Y(S)
      INTEGER           OB(S), RL(S)
      LOGICAL           X(NROWXR,N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A2, SUM, SUM2, SUM3
      INTEGER           I, I2, I3, IFAILA, J, J3, K, L
      LOGICAL           CHANGE, TEMP
C     .. Local Arrays ..
      CHARACTER*80      REC(6)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      EXTERNAL          G01ECF
C     .. External Subroutines ..
      EXTERNAL          G11SAU, G11SAY, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, MIN, DBLE
C     .. Executable Statements ..
C
C     RECODE ITEMS WHICH GIVE NEGATIVE FACTOR LOADINGS
C     AND 'REVERSE' THE SIGN OF THE CORRESPONDING PARAMETERS
C
      CHANGE = .FALSE.
      DO 60 J = 1, N2
         IF (ALPHA(J).LT.0.0D0) THEN
            CHANGE = .TRUE.
            DO 20 L = 1, S
               X(L,J) = .NOT. (X(L,J))
   20       CONTINUE
            ALPHA(J) = -ALPHA(J)
            A(J) = -A(J)
            C(J) = -C(J)
            IF (GPROB) THEN
               GAMMA(J) = -GAMMA(J)
            ELSE
               GAMMA(J) = 1.0D0 - GAMMA(J)
            END IF
C
            G(2*J-1) = -G(2*J-1)
            G(2*J) = -G(2*J)
C
            DO 40 I = 1, N2
               IF (I.NE.J) THEN
                  VAR(2*I-1,2*J-1) = -VAR(2*I-1,2*J-1)
                  VAR(2*I,2*J) = -VAR(2*I,2*J)
                  VAR(2*I-1,2*J) = -VAR(2*I-1,2*J)
                  VAR(2*I,2*J-1) = -VAR(2*I,2*J-1)
               END IF
   40       CONTINUE
C
         END IF
   60 CONTINUE
C
      IF ((CHANGE) .AND. ((ISHOW.EQ.1) .OR. (ISHOW.EQ.4)
     *    .OR. (ISHOW.EQ.5) .OR. (ISHOW.EQ.7))) THEN
         WRITE (REC,FMT=99997)
         CALL X04BAY(LP,2,REC)
      END IF
C
C     CALCULATE THE FACTOR SCORES
C
      CALL G11SAU(Y,XL,X,P,S,N2,AX,XN,Q,ALPHA,GAMMA,A,C,GPROB,PHI,
     *            NROWXR)
C
C     CORRECT THE EXPECTED FREQUENCIES FOR THE S DIFFERENT
C     SCORE PATTERNS
C
      DO 80 L = 1, S
         P(L) = P(L)*DBLE(N)
   80 CONTINUE
C
C     CALCULATE OBSERVED AND EXPECTED FIRST AND SECOND
C     ORDER MARGINS
C
      CALL G11SAY(X,N2,S,RL,LP,N,AX,PHI,Q,ISHOW,OBS,EXPP,NROWXR,IA)
C
C     COMPUTE NUMBER OF POSITIVE RESPONSES FOR EACH SCORE PATTERN
C     AND STORE IN ARRAY OB
C
      DO 120 L = 1, S
         I2 = 0
         DO 100 J = 1, N2
            IF (X(L,J)) I2 = I2 + 1
  100    CONTINUE
         OB(L) = I2
  120 CONTINUE
C
C     DISPLAY OBSERVED AND EXPECTED FREQUENCIES AND FACTOR SCORES
C
C     FIRST SORT THE THETA SCORES INTO ASCENDING ORDER
C
      DO 180 I = 1, S - 1
         DO 160 J = I + 1, S
            IF (Y(J).LT.Y(I)) THEN
C
C              INTERCHANGE ALL SCORES, FREQUENCIES AND SCORE PATTERNS
C
               I2 = OB(I)
               OB(I) = OB(J)
               OB(J) = I2
               A2 = P(I)
               P(I) = P(J)
               P(J) = A2
               J3 = RL(I)
               RL(I) = RL(J)
               RL(J) = J3
C
               A2 = Y(I)
               Y(I) = Y(J)
               Y(J) = A2
               IF ( .NOT. GPROB) THEN
                  A2 = XL(I)
                  XL(I) = XL(J)
                  XL(J) = A2
               END IF
               DO 140 K = 1, N2
                  TEMP = X(I,K)
                  X(I,K) = X(J,K)
                  X(J,K) = TEMP
  140          CONTINUE
C
            END IF
  160    CONTINUE
  180 CONTINUE
C
      IF ((ISHOW.LE.2) .OR. (ISHOW.EQ.4)) GO TO 300
C
      IF ( .NOT. GPROB) THEN
C
         WRITE (REC,FMT=99996)
         CALL X04BAY(LP,6,REC)
C
         DO 220 L = 1, S
            IF (N2.EQ.23) THEN
               WRITE (REC,FMT=99999) RL(L), P(L), Y(L), XL(L), OB(L),
     *           (X(L,J),J=1,N2)
               CALL X04BAF(LP,REC(1))
            ELSE
               WRITE (REC,FMT=99995) RL(L), P(L), Y(L), XL(L), OB(L),
     *           (X(L,J),J=1,MIN(23,N2))
               CALL X04BAF(LP,REC(1))
               DO 200 J = 24, N2, 23
                  WRITE (REC,FMT=99990) (X(L,I),I=J,MIN(J+22,N2))
                  CALL X04BAF(LP,REC(1))
  200          CONTINUE
            END IF
  220    CONTINUE
C
      ELSE
C
         WRITE (REC,FMT=99994)
         CALL X04BAY(LP,6,REC)
C
         DO 260 L = 1, S
            IF (N2.EQ.35) THEN
               WRITE (REC,FMT=99998) RL(L), P(L), Y(L), OB(L),
     *           (X(L,J),J=1,N2)
               CALL X04BAF(LP,REC(1))
            ELSE
               WRITE (REC,FMT=99993) RL(L), P(L), Y(L), OB(L),
     *           (X(L,J),J=1,MIN(35,N2))
               CALL X04BAF(LP,REC(1))
               DO 240 J = 36, N2, 35
                  WRITE (REC,FMT=99989) (X(L,I),I=J,MIN(J+34,N2))
                  CALL X04BAF(LP,REC(1))
  240          CONTINUE
            END IF
  260    CONTINUE
C
      END IF
C
      SUM = 0.0D0
      J3 = 0
      DO 280 L = 1, S
         SUM = SUM + P(L)
         J3 = J3 + RL(L)
  280 CONTINUE
      WRITE (REC,FMT=99992) J3, SUM
      CALL X04BAY(LP,3,REC)
C
  300 IF ( .NOT. CHISQR) RETURN
C
C     CALCULATE LIKELIHOOD RATIO STATISTIC OF OBSERVED AND
C     EXPECTED FREQUENCIES
C
      SUM = 0.0D0
      I = 0
      I2 = 0
      SUM2 = 0.0D0
C
      TEMP = .TRUE.
      DO 340 L = 1, S
         IF (RL(L).EQ.0) GO TO 340
         SUM2 = SUM2 + P(L)
         I2 = I2 + RL(L)
         IF (SUM2.GE.5.0D0) THEN
            IF (L.LT.S) THEN
               SUM3 = 0.0D0
               I3 = 0
               DO 320 J = L + 1, S
                  SUM3 = SUM3 + P(J)
                  I3 = I3 + RL(J)
  320          CONTINUE
               IF (SUM3.LT.5.0D0) THEN
                  SUM2 = SUM2 + SUM3
                  I2 = I2 + I3
                  TEMP = .FALSE.
               END IF
            END IF
C
            SUM = SUM + (DBLE(I2)*LOG(DBLE(I2)/SUM2))
            I = I + 1
            SUM2 = 0.0D0
            I2 = 0
            IF ( .NOT. TEMP) GO TO 360
C
         END IF
C
  340 CONTINUE
  360 I2 = 1
      DO 380 J = 1, N2
         I2 = I2*2
         IF (I2.GT.I) GO TO 400
  380 CONTINUE
      I = I - 1
  400 I = I - 2*N2
      IF ((I.LE.0) .AND. (IERROR.EQ.0)) IERROR = 4
      CHI = 2.0D0*SUM
      IF ((CHI.GT.0.0D0) .AND. (I.GT.0)) THEN
         IFAILA = 1
         SIGLEV = G01ECF('U',CHI,DBLE(I),IFAILA)
      ELSE
         SIGLEV = 0.0D0
      END IF
      IF ((ISHOW.EQ.3) .OR. (ISHOW.GE.5)) THEN
         WRITE (REC,FMT=99991) CHI, SIGLEV, I
         CALL X04BAY(LP,4,REC)
      END IF
C
      IDF = I
      RETURN
C
99999 FORMAT (' ',I6,F13.3,F8.3,F10.3,I8,4X,23L1)
99998 FORMAT (' ',I6,F13.3,F8.3,I6,4X,35L1)
99997 FORMAT (/' THE CODING HAS BEEN REVERSED ON ITEMS WHOSE FACTOR LO',
     *  'ADINGS WERE NEGATIVE')
99996 FORMAT (/' ',/' OBSERVED',3X,'EXPECTED',3X,'THETA',3X,'COMPONENT',
     *  3X,'RAW',4X,'SCORE',/' FREQUENCY',2X,'FREQUENCY',2X,'SCORE',3X,
     *  'SCORE',7X,'SCORE',2X,'PATTERN',/' ---------',2X,'---------',2X,
     *  '-----',3X,'---------',3X,'-----',2X,'-------',/)
99995 FORMAT (' ',I6,F13.3,F8.3,F10.3,I8,4X,23L1)
99994 FORMAT (/' ',/' OBSERVED',3X,'EXPECTED',3X,'THETA',3X,'RAW',4X,
     *  'SCORE',/' FREQUENCY',2X,'FREQUENCY',2X,'SCORE',3X,'SCORE',2X,
     *  'PATTERN',/' ---------',2X,'---------',2X,'-----',3X,'-----',2X,
     *  '-------',/)
99993 FORMAT (' ',I6,F13.3,F8.3,I6,4X,35L1)
99992 FORMAT (' ---------',2X,'---------',/I7,F13.3,/)
99991 FORMAT (/' LIKELIHOOD RATIO GOODNESS OF FIT STATISTIC = ',F10.3,
     *  /25X,'SIGNIFICANCE LEVEL = ',F10.3,/' (BASED ON ',I3,' DEGREES',
     *  ' OF FREEDOM)')
99990 FORMAT (50X,23L1)
99989 FORMAT (38X,35L1)
      END
