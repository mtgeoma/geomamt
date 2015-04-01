      SUBROUTINE G13DCT(LP,IFAIL,LOGL,P,Q,K,X,N4,DISP,IDISP,MEAN,QQ,IK,
     *                  ISHOW,N,V)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     DISPLAY PARAMETER ESTIMATES AND THEIR STANDARD ERRORS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  LOGL
      INTEGER           IDISP, IFAIL, IK, ISHOW, K, LP, N, N4, P, Q
      LOGICAL           MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  DISP(IDISP,N4), QQ(IK,K), V(IK,N), X(N4)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA
      INTEGER           I, I2, II, J, K2, L
C     .. Local Arrays ..
      CHARACTER*120     REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
      WRITE (REC,FMT=99999) IFAIL
      CALL X04BAY(LP,2,REC)
      WRITE (REC,FMT=99998) LOGL
      CALL X04BAY(LP,2,REC)
C
C     PRINT OUT AR PARAMETER ESTIMATES AND THEIR STANDARD ERRORS
C
      IF (P.EQ.0) GO TO 140
C
      IF (K.EQ.1) THEN
         WRITE (REC,FMT=99997)
         CALL X04BAY(LP,3,REC)
      ELSE
         WRITE (REC,FMT=99996)
         CALL X04BAY(LP,3,REC)
      END IF
C
      DO 120 L = 1, P
         DO 100 I = 1, K
            I2 = (L-1)*K*K + (I-1)*K
            IF (K.EQ.1) THEN
               IF (L.LE.9) THEN
                  WRITE (REC,FMT=99995) L, (X(I2+J),J=1,K)
                  CALL X04BAY(LP,2,REC)
               END IF
               IF (L.GT.9) THEN
                  WRITE (REC,FMT=99994) L, (X(I2+J),J=1,K)
                  CALL X04BAY(LP,2,REC)
               END IF
            ELSE
               IF ((I.EQ.1) .AND. (L.LE.9)) THEN
                  WRITE (REC,FMT=99993) L, (X(I2+J),J=1,MIN(10,K))
                  CALL X04BAY(LP,2,REC)
                  DO 20 J = 11, K, 10
                     WRITE (REC,FMT=99991) (X(I2+II),II=J,MIN(J+9,K))
                     CALL X04BAY(LP,2,REC)
   20             CONTINUE
               END IF
               IF ((I.EQ.1) .AND. (L.GT.9)) THEN
                  WRITE (REC,FMT=99992) L, (X(I2+J),J=1,MIN(10,K))
                  CALL X04BAY(LP,2,REC)
                  DO 40 J = 11, K, 10
                     WRITE (REC,FMT=99991) (X(I2+II),II=J,MIN(J+9,K))
                     CALL X04BAY(LP,2,REC)
   40             CONTINUE
               END IF
               IF (I.GT.1) THEN
                  DO 60 J = 1, K, 10
                     WRITE (REC,FMT=99991) (X(I2+II),II=J,MIN(J+9,K))
                     CALL X04BAY(LP,2,REC)
   60             CONTINUE
               END IF
            END IF
            IF (K.EQ.1) THEN
               WRITE (REC,FMT=99989) (DISP(I2+J,I2+J),J=1,K)
               CALL X04BAF(LP,REC(1))
            END IF
            IF (K.GT.1) THEN
               WRITE (REC,FMT=99987) (DISP(I2+J,I2+J),J=1,MIN(10,K))
               CALL X04BAF(LP,REC(1))
               DO 80 J = 11, K, 10
                  WRITE (REC,FMT=99987) (DISP(I2+II,I2+II),II=J,MIN(J+9,
     *              K))
                  CALL X04BAF(LP,REC(1))
   80          CONTINUE
            END IF
  100    CONTINUE
  120 CONTINUE
C
C     NOW PRINT OUT MA PARAMETER ESTIMATES AND THEIR STANDARD ERRORS
C
  140 IF (Q.EQ.0) GO TO 280
C
      IF (K.EQ.1) THEN
         WRITE (REC,FMT=99986)
         CALL X04BAY(LP,3,REC)
      ELSE
         WRITE (REC,FMT=99985)
         CALL X04BAY(LP,3,REC)
      END IF
C
      DO 260 L = 1, Q
         DO 240 I = 1, K
            I2 = P*K*K + (L-1)*K*K + (I-1)*K
            IF (K.EQ.1) THEN
               IF (L.LE.9) THEN
                  WRITE (REC,FMT=99984) L, (X(I2+J),J=1,K)
                  CALL X04BAY(LP,2,REC)
               END IF
               IF (L.GT.9) THEN
                  WRITE (REC,FMT=99983) L, (X(I2+J),J=1,K)
                  CALL X04BAY(LP,2,REC)
               END IF
            ELSE
               IF ((I.EQ.1) .AND. (L.LE.9)) THEN
                  WRITE (REC,FMT=99982) L, (X(I2+J),J=1,MIN(10,K))
                  CALL X04BAY(LP,2,REC)
                  DO 160 J = 11, K, 10
                     WRITE (REC,FMT=99991) (X(I2+II),II=J,MIN(J+9,K))
                     CALL X04BAY(LP,2,REC)
  160             CONTINUE
               END IF
               IF ((I.EQ.1) .AND. (L.GT.9)) THEN
                  WRITE (REC,FMT=99981) L, (X(I2+J),J=1,MIN(10,K))
                  CALL X04BAY(LP,2,REC)
                  DO 180 J = 11, K, 10
                     WRITE (REC,FMT=99991) (X(I2+II),II=J,MIN(J+9,K))
                     CALL X04BAY(LP,2,REC)
  180             CONTINUE
               END IF
               IF (I.GT.1) THEN
                  DO 200 J = 1, K, 10
                     WRITE (REC,FMT=99991) (X(I2+II),II=1,MIN(J+9,K))
                     CALL X04BAY(LP,2,REC)
  200             CONTINUE
               END IF
            END IF
            IF (K.EQ.1) THEN
               WRITE (REC,FMT=99989) (DISP(I2+J,I2+J),J=1,K)
               CALL X04BAF(LP,REC(1))
            END IF
            IF (K.GT.1) THEN
               DO 220 J = 1, K, 10
                  WRITE (REC,FMT=99987) (DISP(I2+II,I2+II),II=1,MIN(J+9,
     *              K))
                  CALL X04BAF(LP,REC(1))
  220          CONTINUE
            END IF
  240    CONTINUE
  260 CONTINUE
C
  280 IF ( .NOT. MEAN) GO TO 340
      WRITE (REC,FMT=99980)
      CALL X04BAY(LP,3,REC)
      J = (P+Q)*K*K
      AA = DISP(J+1,J+1)
      DO 300 I = 2, K
         AA = MAX(DISP(J+I,J+I),AA)
  300 CONTINUE
C
      IF (K.EQ.1) THEN
         WRITE (REC,FMT=99979) (X((P+Q)*K*K+I),I=1,K)
         CALL X04BAY(LP,2,REC)
         IF (AA.LE.9.999D0) THEN
            WRITE (REC,FMT=99989) (DISP((P+Q)*K*K+I,(P+Q)*K*K+I),I=1,K)
         ELSE
            WRITE (REC,FMT=99990) (DISP(J+I,J+I),I=1,K)
         END IF
         CALL X04BAF(LP,REC(1))
      ELSE
         DO 320 I = 1, K, 10
            WRITE (REC,FMT=99991) (X((P+Q)*K*K+J),J=I,MIN(I+9,K))
            CALL X04BAY(LP,2,REC)
            IF (AA.LE.9.999D0) THEN
               WRITE (REC,FMT=99987) (DISP((P+Q)*K*K+J,(P+Q)*K*K+J),J=I,
     *           MIN(I+9,K))
            ELSE
               WRITE (REC,FMT=99988) (DISP((P+Q)*K*K+J,(P+Q)*K*K+J),J=I,
     *           MIN(I+9,K))
            END IF
            CALL X04BAF(LP,REC(1))
  320    CONTINUE
      END IF
C
C     PRINT OUT SIGMA MATRIX
C
  340 IF (K.EQ.1) THEN
         WRITE (REC,FMT=99978)
         CALL X04BAY(LP,3,REC)
      ELSE
         WRITE (REC,FMT=99977)
         CALL X04BAY(LP,3,REC)
C
      END IF
C
      DO 380 I = 1, K
         IF (K.EQ.1) THEN
            WRITE (REC,FMT=99979) QQ(1,1)
            CALL X04BAY(LP,2,REC)
         ELSE
            DO 360 J = 1, I, 10
               WRITE (REC,FMT=99991) (QQ(I,II),II=J,MIN(J+9,I))
               CALL X04BAY(LP,2,REC)
  360       CONTINUE
         END IF
  380 CONTINUE
C
      IF (ISHOW.EQ.1) RETURN
C
C     PRINT OUT RESIDUAL SERIES
C
      DO 460 I = 1, K
         WRITE (REC,FMT=99976) I
         CALL X04BAY(LP,3,REC)
         DO 400 J = 1, N/8
            K2 = (J-1)*8
            WRITE (REC,FMT=99975) (K2+I2,I2=1,8)
            CALL X04BAY(LP,2,REC)
            WRITE (REC,FMT=99974) (V(I,K2+I2),I2=1,8)
            CALL X04BAF(LP,REC(1))
  400    CONTINUE
         J = 8*(N/8)
         IF (N.GT.J) THEN
            DO 420 II = J + 1, N, 8
               WRITE (REC,FMT=99975) (I2,I2=II,MIN(II+7,N))
               CALL X04BAY(LP,2,REC)
  420       CONTINUE
         END IF
         IF (N.GT.J) THEN
            DO 440 II = J + 1, N, 8
               WRITE (REC,FMT=99974) (V(I,I2),I2=II,MIN(II+7,N))
               CALL X04BAF(LP,REC(1))
  440       CONTINUE
         END IF
  460 CONTINUE
C
      RETURN
C
99999 FORMAT (/' VALUE OF IFAIL PARAMETER ON EXIT FROM G13DCF = ',I3)
99998 FORMAT (/' VALUE OF LOG LIKELIHOOD FUNCTION ON EXIT = ',D12.5)
99997 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATES OF AR PARAMETERS',/' ---',
     *  '------------------------------------------')
99996 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATES OF AR PARAMETER MATRICES',
     *  /' -----------------------------------------------------')
99995 FORMAT (/' PHI(',I1,')    : ',F8.3)
99994 FORMAT (/' PHI(',I2,')   : ',F8.3)
99993 FORMAT (/' PHI(',I1,')    ROW-WISE : ',10F8.3)
99992 FORMAT (/' PHI(',I2,')   ROW-WISE : ',10F8.3)
99991 FORMAT (/22X,10F8.3)
99990 FORMAT (14X,10('(',F6.3,')',:))
99989 FORMAT (14X,10(' (',F5.3,')',:))
99988 FORMAT (23X,10('(',F6.3,')',:))
99987 FORMAT (23X,10(' (',F5.3,')',:))
99986 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATES OF MA PARAMETERS',/' ---',
     *  '------------------------------------------')
99985 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATES OF MA PARAMETER MATRICES',
     *  /' -----------------------------------------------------')
99984 FORMAT (/' THETA(',I1,')  : ',F8.3)
99983 FORMAT (/' THETA(',I2,') : ',F8.3)
99982 FORMAT (/' THETA(',I1,')  ROW-WISE : ',10F8.3)
99981 FORMAT (/' THETA(',I2,') ROW-WISE : ',10F8.3)
99980 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATE OF PROCESS MEAN',/' -----',
     *  '--------------------------------------')
99979 FORMAT (/13X,F8.3)
99978 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATE OF SIGMA',/' ------------',
     *  '------------------------')
99977 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATE OF SIGMA MATRIX',/' -----',
     *  '--------------------------------------')
99976 FORMAT (/11X,'RESIDUAL SERIES NUMBER ',I2,/11X,'----------------',
     *  '---------')
99975 FORMAT (/'   T ',I5,7I7)
99974 FORMAT (' V(T)',8F7.2)
      END
