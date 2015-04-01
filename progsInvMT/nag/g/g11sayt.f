      SUBROUTINE G11SAY(X,N2,S,RL,LP,N,AX,PHI,Q,ISHOW,OBS,EXPP,NROWXR,
     *                  IA)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATE OBSERVED AND EXPECTED FIRST AND SECOND
C     ORDER MARGINS
C
C     .. Scalar Arguments ..
      INTEGER           IA, ISHOW, LP, N, N2, NROWXR, Q, S
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(20), EXPP(IA,N2), OBS(IA,N2), PHI(N2,20)
      INTEGER           RL(S)
      LOGICAL           X(NROWXR,N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, ISUM, IT, ITIMES, J, K, L
      CHARACTER*72      ST
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, MOD, DBLE
C     .. Executable Statements ..
C
C     CALCULATE OBSERVED NUMBER OF INDIVIDUALS WHO RESPOND CORRECTLY
C     TO INDIVIDUAL AND PAIRS OF ITEMS
C
      DO 80 I = 1, N2
         DO 60 J = 1, I
            IF (I.EQ.J) THEN
C
C              COUNT NUMBER OF INDIVIDUALS ANSWERING ITEM I
C              CORRECTLY
C
               ISUM = 0
               DO 20 L = 1, S
                  IF (X(L,I)) ISUM = ISUM + RL(L)
   20          CONTINUE
               OBS(I,I) = (100.0D0/DBLE(N))*ISUM
            ELSE
C
C              COUNT NUMBER OF INDIVIDUALS ANSWERING ITEMS I AND J
C              CORRECTLY
C
               ISUM = 0
               DO 40 L = 1, S
                  IF (X(L,I) .AND. X(L,J)) ISUM = ISUM + RL(L)
   40          CONTINUE
               OBS(I,J) = (100.0D0/DBLE(N))*ISUM
C
            END IF
C
   60    CONTINUE
   80 CONTINUE
C
C     EVALUATE 1ST AND 2 ND ORDER MARGINS
C
      DO 160 I = 1, N2
         IF (I.NE.1) THEN
            DO 120 J = 1, I - 1
               SUM = 0.0D0
               DO 100 K = 1, Q
                  SUM = SUM + PHI(I,K)*PHI(J,K)*AX(K)
  100          CONTINUE
               EXPP(I,J) = SUM*100.0D0
  120       CONTINUE
         END IF
         SUM = 0.0D0
         DO 140 K = 1, Q
            SUM = SUM + PHI(I,K)*AX(K)
  140    CONTINUE
         EXPP(I,I) = SUM*100.0D0
  160 CONTINUE
C
      IF ((ISHOW.LE.1) .OR. (ISHOW.EQ.3) .OR. (ISHOW.EQ.5)) RETURN
C
C     DISPLAY FIRST AND SECOND ORDER OBSERVED AND EXPECTED VALUES
C
      WRITE (REC,FMT=99999)
      CALL X04BAY(LP,5,REC)
C
      ITIMES = N2/12
      IF (MOD(N2,12).NE.0) ITIMES = ITIMES + 1
C
      DO 220 J = 1, ITIMES
C
         DO 180 I = 1, 72
            ST(I:I) = ' '
            IF (MOD(I,6).EQ.0) THEN
               ST(I:I) = '-'
               IF ((J.GT.1) .OR. (I.GT.54)) ST(I-1:I-1) = '-'
            END IF
  180    CONTINUE
C
         IT = MIN(N2,12*(J-1)+12)
         WRITE (REC,FMT=99998) (I,I=12*(J-1)+1,MIN(J*12,IT))
         CALL X04BAY(LP,5,REC)
         I = (IT-12*(J-1))*6
         WRITE (REC,FMT=99997) ST(1:I)
         CALL X04BAF(LP,REC(1))
C
         DO 200 I = 12*(J-1) + 1, N2
            WRITE (REC,FMT=99996) I, (EXPP(I,K),K=12*(J-1)+1,MIN(I,12*J)
     *        )
            CALL X04BAF(LP,REC(1))
            WRITE (REC,FMT=99995) (OBS(I,K),K=12*(J-1)+1,MIN(I,12*J))
            CALL X04BAF(LP,REC(1))
  200    CONTINUE
C
  220 CONTINUE
C
      RETURN
C
99999 FORMAT (/' ',/' EXPECTED (AND OBSERVED) PERCENTAGE OF CASES PROD',
     *  'UCING',/' POSITIVE RESPONSES FOR INDIVIDUAL AND PAIRS OF ITEMS'
     *  ,/' -----------------------------------------------------')
99998 FORMAT (/41X,'ITEM',/41X,'----',//' ITEM ',12I6)
99997 FORMAT (' ---- ',A)
99996 FORMAT (' ',I3,3X,12F6.1)
99995 FORMAT (8X,12('(',F4.1,')',:))
      END
