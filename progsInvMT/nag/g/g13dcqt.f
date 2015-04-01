      SUBROUTINE G13DCQ(N2,X,F,G,ISTATE,GPJNRM,COND,POSDEF,NITER,NF,IW,
     *                  LIW,W2,LW)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     AUTOMATIC MONITORING ROUTINE (MONIT) FOR E04JBL
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, F, GPJNRM
      INTEGER           LIW, LW, N2, NF, NITER
      LOGICAL           POSDEF
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N2), W2(LW), X(N2)
      INTEGER           ISTATE(N2), IW(LIW)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ADDLOG, CONDD, LMAX, NORM
      INTEGER           ITN, K, K4, K5, K6, KR, LEW6, LEW7, LP, LW1,
     *                  LW11, LW12, LW13, LW14, LW16, LW17, LW18, LW2,
     *                  LW3, LW4, LW5, LW8, P, Q
      LOGICAL           FULLP, FULLQ, MEAN, NOPRIN
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, I3, I4, I9, J, L, N
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          G13DCV, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /BG13DC/ADDLOG, LMAX, CONDD, NORM, K, P, Q, K4,
     *                  K5, LW14, LW18, LP, K6, ITN, LEW6, LEW7, LW2,
     *                  LW1, LW3, LW4, LW5, LW8, LW11, LW12, LW13, LW16,
     *                  LW17, KR, MEAN, NOPRIN, FULLP, FULLQ
C     .. Executable Statements ..
      CONDD = COND
      NORM = GPJNRM
      ITN = NITER
      IF (NOPRIN) RETURN
      WRITE (REC,FMT=99999) NITER, COND
      CALL X04BAY(LP,4,REC)
      WRITE (REC,FMT=99998) NF
      CALL X04BAY(LP,2,REC)
      WRITE (REC,FMT=99997) - F - ADDLOG
      CALL X04BAY(LP,3,REC)
C
C     RECONSTRUCT SIGMA MATRIX
C
      DO 40 I = 1, K
         DO 20 J = 1, I
            SUM = 0.0D0
            IF (I.EQ.J) SUM = 1.0D0
            W2(LW14-1+(I-1)*K+J) = SUM
            W2(LW14-1+(J-1)*K+I) = SUM
C
   20    CONTINUE
   40 CONTINUE
C
      IF (P.EQ.0) GO TO 200
C
      IF (FULLP) THEN
C
C        RECOVER THE ORIGINAL AR PARAMETERS FROM THE
C        TRANSFORMED VARIABLES
C
         I9 = 0
         CALL G13DCV(X,P,W2(LW1),K,N2,W2(LW14),W2(LEW6),W2(LW3),W2(LW4),
     *               W2(LW5),W2(LW17),W2(LW8),W2(LW11),W2(LW16),W2(LW12)
     *               ,W2(LW13),W2(LEW7),KR,I9,I4)
         DO 100 L = 1, P
            DO 80 J = 1, K
               DO 60 I = 1, K
                  W2(LEW6-1+(L-1)*K*K+(I-1)*K+J) = W2(LW1-1+(L-1)
     *              *K*K+(J-1)*K+I)
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
C
      ELSE
C
         DO 160 L = 1, P
            DO 140 I = 1, K
               DO 120 J = 1, K
                  W2(LEW6-1+(L-1)*K*K+(I-1)*K+J) = X((L-1)*K*K+(I-1)
     *              *K+J)
  120          CONTINUE
  140       CONTINUE
  160    CONTINUE
C
      END IF
C
      IF ((K4.EQ.5) .AND. (P.GT.0)) THEN
         WRITE (REC,FMT=99996) (W2(LEW6-1+I),I=1,5)
         CALL X04BAF(LP,REC(1))
      END IF
      IF ((K4.NE.5) .AND. (P.GT.0)) THEN
         WRITE (REC,FMT=99995) (W2(LEW6-1+I),I=1,MIN(5,K4))
         CALL X04BAF(LP,REC(1))
         DO 180 I = 6, K4, 5
            WRITE (REC,FMT=99985) (W2(LEW6-1+J),J=I,MIN(I+4,K4))
            CALL X04BAF(LP,REC(1))
  180    CONTINUE
      END IF
C
  200 IF (Q.EQ.0) GO TO 380
C
      IF (FULLQ) THEN
C
C        RECOVER THE ORIGINAL MA PARAMETERS FROM THE TRANSFORMED
C        VARIABLES
C
         CALL G13DCV(X,Q,W2(LW2),K,N2,W2(LW14),W2(LEW6),W2(LW3),W2(LW4),
     *               W2(LW5),W2(LW17),W2(LW8),W2(LW11),W2(LW16),W2(LW12)
     *               ,W2(LW13),W2(LEW7),KR,K4,I4)
         DO 260 L = 1, Q
            DO 240 J = 1, K
               DO 220 I = 1, K
                  W2(LEW6-1+(L-1)*K*K+(I-1)*K+J) = W2(LW2-1+(L-1)
     *              *K*K+(J-1)*K+I)
  220          CONTINUE
  240       CONTINUE
  260    CONTINUE
C
      ELSE
C
         DO 320 L = 1, Q
            DO 300 I = 1, K
               DO 280 J = 1, K
                  W2(LEW6-1+(L-1)*K*K+(I-1)*K+J) = X(P*K*K+(L-1)
     *              *K*K+(I-1)*K+J)
  280          CONTINUE
  300       CONTINUE
  320    CONTINUE
C
      END IF
C
      IF ((K5-K4.EQ.5) .AND. (Q.GT.0)) THEN
         WRITE (REC,FMT=99994) (W2(LEW6-1+I),I=1,MIN(5,Q*K*K))
         CALL X04BAF(LP,REC(1))
         DO 340 I = 6, Q*K*K, 5
            WRITE (REC,FMT=99985) (W2(LEW6-1+J),J=I,MIN(I+4,Q*K*K))
            CALL X04BAF(LP,REC(1))
  340    CONTINUE
      END IF
      IF ((K5-K4.NE.5) .AND. (Q.GT.0)) THEN
         WRITE (REC,FMT=99993) (W2(LEW6-1+I),I=1,MIN(5,Q*K*K))
         CALL X04BAF(LP,REC(1))
         DO 360 I = 6, Q*K*K, 5
            WRITE (REC,FMT=99985) (W2(LEW6-1+J),J=I,MIN(I+4,Q*K*K))
            CALL X04BAF(LP,REC(1))
  360    CONTINUE
      END IF
  380 IF (MEAN .AND. (K.EQ.5)) THEN
         WRITE (REC,FMT=99992) (X(K5+I)+W2(LW18-1+I),I=1,K)
         CALL X04BAF(LP,REC(1))
      END IF
      IF (MEAN .AND. (K.NE.5)) THEN
         WRITE (REC,FMT=99991) (X(K5+I)+W2(LW18-1+I),I=1,MIN(5,K))
         CALL X04BAF(LP,REC(1))
         DO 400 I = 6, K, 5
            WRITE (REC,FMT=99985) (X(K5+J)+W2(LW18-1+J),J=I,MIN(I+4,K))
            CALL X04BAF(LP,REC(1))
  400    CONTINUE
      END IF
C
C     RECONSTRUCT SIGMA MATRIX
C
      DO 460 I = 1, K
         DO 440 J = 1, I
            SUM = 0.0D0
            DO 420 I3 = 1, K
               IF (MIN(I,J).GE.I3) SUM = SUM + X(K6+(I-1)*I/2+I3)
     *                                   *X(K6+(J-1)*J/2+I3)
  420       CONTINUE
            W2(LW14-1+(I-1)*K+J) = LMAX*LMAX*SUM
            W2(LW14-1+(J-1)*K+I) = LMAX*LMAX*SUM
C
  440    CONTINUE
  460 CONTINUE
C
      DO 500 I = 1, K
         IF (I.EQ.1) THEN
            WRITE (REC,FMT=99990) W2(LW14)
            CALL X04BAF(LP,REC(1))
         END IF
         IF ((I.GT.1) .AND. (I.NE.5)) THEN
            WRITE (REC,FMT=99989) (W2(LW14-1+(I-1)*K+J),J=1,MIN(5,I))
            CALL X04BAF(LP,REC(1))
            DO 480 J = 6, I, 5
               WRITE (REC,FMT=99985) (W2(LW14-1+(I-1)*K+N),N=J,MIN(J+4,
     *           I))
               CALL X04BAF(LP,REC(1))
  480       CONTINUE
         END IF
         IF (I.EQ.5) THEN
            WRITE (REC,FMT=99988) (W2(LW14-1+(I-1)*K+J),J=1,I)
            CALL X04BAF(LP,REC(1))
         END IF
  500 CONTINUE
C
      WRITE (REC,FMT=99987) GPJNRM
      CALL X04BAY(LP,2,REC)
      WRITE (REC,FMT=99986)
      CALL X04BAY(LP,2,REC)
C
      RETURN
C
99999 FORMAT (/' ITERATION NUMBER = ',I5,//' ESTIMATED CONDITION NUMBE',
     *  'R OF HESSIAN MATRIX = ',D12.3)
99998 FORMAT (/' NUMBER OF LIKELIHOOD EVALUATIONS MADE SO FAR = ',I8)
99997 FORMAT (/' VALUE OF LOG LIKELIHOOD FUNCTION = ',D12.5,/)
99996 FORMAT (' AR PARAMETERS : ',5D11.3)
99995 FORMAT (' AR PARAMETERS : ',5D11.3)
99994 FORMAT (' MA PARAMETERS : ',5D11.3)
99993 FORMAT (' MA PARAMETERS : ',5D11.3)
99992 FORMAT (' MEAN VECTOR   : ',5D11.3)
99991 FORMAT (' MEAN VECTOR   : ',5D11.3)
99990 FORMAT (' SIGMA MATRIX  : ',D11.3)
99989 FORMAT (17X,5D11.3)
99988 FORMAT (17X,5D11.3)
99987 FORMAT (/' NORM OF GRADIENT VECTOR = ',D12.5)
99986 FORMAT (/' *****************************************************',
     *  '********************')
99985 FORMAT (17X,5D11.3)
      END
