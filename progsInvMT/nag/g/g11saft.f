      SUBROUTINE G11SAF(N2,N,GPROB,S,X,NROWXR,RL,A,C,IPRINT,CGETOL,
     *                  MAXIT,CHISQR,ISHOW,NITER,ALPHA,GAMMA,VAR,IAA,G,
     *                  EXPP,IA,OBS,P,Y,XL,OB,LL,CHI,IDF,SIGLEV,W,LW,
     *                  IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13A REVISED. IER-634 (APR 1988).
C
C     MAXIMUM LIKELIHOOD ESTIMATION OF ITEM PARAMETERS VIA THE
C     E-M ALGORITHM FOR THE FOLLOWING MODELS
C
C          (1) LOGIT
C          (2) PROBIT
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G11SAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CGETOL, CHI, LL, SIGLEV
      INTEGER           IA, IAA, IDF, IFAIL, IPRINT, ISHOW, LW, MAXIT,
     *                  N, N2, NITER, NROWXR, S
      LOGICAL           CHISQR, GPROB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N2), ALPHA(N2), C(N2), EXPP(IA,N2), G(2*N2),
     *                  GAMMA(N2), OBS(IA,N2), P(S), VAR(IAA,2*N2),
     *                  W(LW), XL(S), Y(S)
      INTEGER           OB(S), RL(S)
      LOGICAL           X(NROWXR,N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A3, GRATOL, H4, LIMIT, ROOTPI, SUM
      INTEGER           I, I2, IAN, IERROR, IEXIT, ISUM, J, K, L, LP,
     *                  LW1, LW2, LW3, LW4, LW5, N3, Q
      LOGICAL           NOMON
C     .. Local Arrays ..
      DOUBLE PRECISION  AX(20), XN(20)
      CHARACTER*80      P01REC(6), REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G11SAT, G11SAW, G11SAX, G11SAZ, X04ABF, X04BAF,
     *                  X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, MIN, SQRT
C     .. Executable Statements ..
C
      CALL X04ABF(0,LP)
C
      IF (CGETOL.LT.SQRT(X02AJF())) THEN
         GRATOL = SQRT(X02AJF())
         WRITE (REC,FMT=99999) CGETOL, GRATOL
         CALL X04BAY(LP,2,REC)
      ELSE
         GRATOL = CGETOL
      END IF
      LIMIT = LOG(X02AMF())
C
C     PARTITION WORK SPACE ARRAY
      N3 = N2 + N2
      LW1 = 1
      LW2 = N2*20 + LW1
      LW3 = LW2 + N2*20
      LW4 = LW3 + N3
      LW5 = LW4 + N2*20
      IAN = N2 + N2 + 1
C
C     TEST FOR ERRORS IN INPUT DATA
C
      IERROR = 1
      IEXIT = IFAIL
      IFAIL = 0
      IF (LW.LT.(4*N2*(N2+16))) THEN
         WRITE (P01REC,FMT=99997) LW, 4*N2*(N2+16)
         IFAIL = P01ABF(IEXIT,IERROR,SRNAME,1,P01REC)
         RETURN
      END IF
      IF (IA.LT.N2) GO TO 520
      IF (IAA.LT.N2+N2) GO TO 520
      IF (MAXIT.LT.1) GO TO 520
C
      I2 = 1
      DO 20 I = 1, N2
         I2 = I2*2
         IF (I2.GT.N) THEN
            IF (S.GT.N) THEN
               GO TO 520
            ELSE
               GO TO 40
            END IF
         END IF
   20 CONTINUE
      IF (S.GT.I2) GO TO 520
C
   40 IF (NROWXR.LT.S) GO TO 520
C
C     INITIALISE OUTPUT CHANNEL NUMBER
C
      NOMON = IPRINT .LT. 0
C
      IF ((ISHOW.LT.0) .OR. (ISHOW.GT.7)) GO TO 520
C
      IF (N.LT.7) GO TO 520
C
      IF (N2.LT.3) GO TO 520
C
      Q = 10
C
      IF (S.LE.2*N2) GO TO 520
      I2 = 0
      DO 60 I = 1, S
         IF (RL(I).LT.0) THEN
            WRITE (P01REC,FMT=99996) I, RL(I)
            IFAIL = P01ABF(IEXIT,IERROR,SRNAME,1,P01REC)
            RETURN
         END IF
         I2 = I2 + RL(I)
   60 CONTINUE
      IF (I2.NE.N) THEN
         WRITE (P01REC,FMT=99995) I2, N
         IFAIL = P01ABF(IEXIT,IERROR,SRNAME,1,P01REC)
         RETURN
      END IF
C
C     TEST WHETHER ANY ROWS OF X ARE REPEATED
C
      DO 120 I = 1, S - 1
         DO 100 J = I + 1, S
            DO 80 K = 1, N2
               IF (X(I,K) .NEQV. X(J,K)) GO TO 100
   80       CONTINUE
            WRITE (P01REC,FMT=99998) I, J
            IFAIL = P01ABF(IEXIT,IERROR,SRNAME,1,P01REC)
            RETURN
  100    CONTINUE
  120 CONTINUE
C
C     IF WE HAVE GOT TO THIS POINT THEN THE INPUT DATA IS
C     FREE FROM ERROR AND WE DO NOT NEED TO RESET IFAIL
C
C
C     TEST WHETHER EVERY ITEM HAS RESPONSES AT BOTH LEVELS
C
      IERROR = 2
      DO 160 J = 1, N2
         ISUM = 0
         DO 140 L = 1, S
            IF (X(L,J)) ISUM = ISUM + RL(L)
  140    CONTINUE
         IF ((ISUM.EQ.0) .OR. (ISUM.EQ.N)) GO TO 520
  160 CONTINUE
C
C     CALCULATE THE ALPHA'S AND GAMMA'S
C
      IF (GPROB) THEN
         DO 180 J = 1, N2
            H4 = 1.0D0/((1.0D0+A(J)*A(J))**0.5D0)
            ALPHA(J) = A(J)*H4
            GAMMA(J) = -C(J)*H4
  180    CONTINUE
      ELSE
         DO 240 J = 1, N2
            ALPHA(J) = A(J)
            IF (C(J).LE.0.0D0) THEN
               IF ((C(J)).LT.(LIMIT)) THEN
                  A3 = 0
                  GO TO 200
               END IF
               A3 = EXP(C(J))
  200          GAMMA(J) = A3/(1.0D0+A3)
            ELSE
               IF ((-C(J)).LT.(LIMIT)) THEN
                  GAMMA(J) = 1.0D0
                  GO TO 220
               END IF
               GAMMA(J) = 1.0D0/(1.0D0+EXP(-C(J)))
  220       END IF
  240    CONTINUE
      END IF
C
C     PRINT OUT ITERATION NUMBER ZERO IF IPRINT .GT. 0
C
      IF (IPRINT.GT.0) THEN
         WRITE (REC,FMT=99990)
         CALL X04BAY(LP,4,REC)
         WRITE (REC,FMT=99989)
         CALL X04BAY(LP,4,REC)
         DO 260 J = 1, N2, 5
            WRITE (REC,FMT=99988) (A(I),I=J,MIN(J+4,N2))
            CALL X04BAF(LP,REC(1))
  260    CONTINUE
         WRITE (REC,FMT=99987)
         CALL X04BAY(LP,4,REC)
         DO 280 J = 1, N2, 5
            WRITE (REC,FMT=99988) (C(I),I=J,MIN(J+4,N2))
            CALL X04BAF(LP,REC(1))
  280    CONTINUE
      END IF
C
C     CALCULATE ONE OVER ROOT TWO PI
C
      ROOTPI = 1.0D0/(SQRT(2.0D0*X01AAF(0.0D0)))
C
C     CALCULATE THE QUADRATURE POINTS AND WEIGHTS
C
      CALL G11SAW(Q,XN,AX,LP,IPRINT,ROOTPI)
C
C     OBTAIN MAXIMUM LIKELIHOOD ESTIMATES OF ITEM PARAMETERS
C
      CALL G11SAX(MAXIT,IPRINT,A,C,XN,AX,Q,N2,S,X,RL,LP,GPROB,GRATOL,
     *            NOMON,ROOTPI,N3,IERROR,W(LW1),G,P,W(LW2),NROWXR,NITER)
C
C     CALCULATE THE ALPHA'S AND GAMMA'S
C
      IF (GPROB) THEN
         DO 300 J = 1, N2
            H4 = 1.0D0/((1.0D0+A(J)*A(J))**0.5D0)
            ALPHA(J) = A(J)*H4
            GAMMA(J) = -C(J)*H4
  300    CONTINUE
      ELSE
         DO 360 J = 1, N2
            ALPHA(J) = A(J)
            IF (C(J).LE.0.0D0) THEN
               IF (C(J).LT.LIMIT) THEN
                  SUM = 0
                  GO TO 320
               END IF
               SUM = EXP(C(J))
  320          GAMMA(J) = SUM/(1.0D0+SUM)
            ELSE
               IF ((-C(J)).LT.(LIMIT)) THEN
                  GAMMA(J) = 1.0D0
                  GO TO 340
               END IF
               GAMMA(J) = 1.0D0/(1.0D0+EXP(-C(J)))
  340       END IF
  360    CONTINUE
      END IF
C
      CALL G11SAZ(A,C,GPROB,N2,LL,VAR,G,XN,AX,Q,S,X,RL,ALPHA,GAMMA,
     *            W(LW1),ROOTPI,N3,W(LW2),P,NROWXR,IAA,W(LW3),W(LW4),
     *            W(LW5),IAN,IERROR)
      IF (IERROR.GT.5) THEN
         DO 400 I = 1, N3
            DO 380 J = 1, N3
               VAR(I,J) = 0.0D0
  380       CONTINUE
  400    CONTINUE
      END IF
C
      IF ((ISHOW.EQ.1) .OR. (ISHOW.EQ.4) .OR. (ISHOW.EQ.5)
     *    .OR. (ISHOW.EQ.7)) THEN
C
         WRITE (REC,FMT=99986) LL
         CALL X04BAY(LP,2,REC)
C
         WRITE (REC,FMT=99985)
         CALL X04BAY(LP,3,REC)
C
         IF (GPROB) THEN
            WRITE (REC,FMT=99984)
            CALL X04BAY(LP,3,REC)
C
            DO 420 I = 1, N2
               WRITE (REC,FMT=99983) I, A(I), ALPHA(I),
     *           SQRT(VAR(2*I-1,2*I-1)), C(I), GAMMA(I),
     *           SQRT(VAR(2*I,2*I))
               CALL X04BAY(LP,2,REC)
  420       CONTINUE
C
         ELSE
C
            WRITE (REC,FMT=99982)
            CALL X04BAY(LP,3,REC)
C
            DO 440 I = 1, N2
               WRITE (REC,FMT=99981) I, ALPHA(I),
     *           SQRT(VAR(2*I-1,2*I-1)), C(I), GAMMA(I),
     *           SQRT(VAR(2*I,2*I))
               CALL X04BAY(LP,2,REC)
  440       CONTINUE
         END IF
C
      END IF
C
C     CORRECT COVARIANCE MATRIX TO A STANDARD ERROR - CORRELATION MATRIX
C
      IF (IERROR.LT.6) THEN
         DO 460 I = 1, N3
            VAR(I,I) = SQRT(VAR(I,I))
  460    CONTINUE
         DO 500 I = 2, N3
            DO 480 J = 1, I - 1
               VAR(I,J) = VAR(I,J)/(VAR(I,I)*VAR(J,J))
               VAR(J,I) = VAR(I,J)
  480       CONTINUE
  500    CONTINUE
      END IF
C
C     PERFORM DIAGNOSTICS
C
      CALL G11SAT(S,N,LP,RL,X,N2,AX,XN,Q,ALPHA,GAMMA,A,C,GPROB,W(LW1),
     *            CHISQR,ISHOW,IERROR,OB,Y,XL,IA,W(LW2),P,NROWXR,OBS,
     *            EXPP,CHI,IDF,VAR,G,IAA,SIGLEV)
C
C     SWITCH IERRORS AROUND
C
  520 IF (IERROR.EQ.4) THEN
         IERROR = 7
         GO TO 540
      END IF
      IF (IERROR.EQ.5) THEN
         IERROR = 4
         GO TO 540
      END IF
      IF (IERROR.EQ.6) THEN
         IERROR = 5
         GO TO 540
      END IF
      IF (IERROR.EQ.7) THEN
         IERROR = 6
         GO TO 540
      END IF
C
  540 IF (ISHOW.GT.0) THEN
         WRITE (REC,FMT=99991) IERROR
         CALL X04BAY(LP,2,REC)
      END IF
C
      IF (IERROR.NE.0) THEN
C
         IF (IERROR.EQ.1) THEN
            WRITE (P01REC,FMT=99994) N2, N, S, ISHOW, IA, IAA, NROWXR,
     *        MAXIT
            IFAIL = P01ABF(IEXIT,IERROR,SRNAME,4,P01REC)
C
         ELSE
C
            IF (IERROR.EQ.2) WRITE (P01REC,FMT=99993)
     *       '** FOR AT LEAST ONE OF THE IP ITEMS THE RESPONSES ARE ALL'
     *          , ' AT THE SAME LEVEL'
            IF (IERROR.EQ.3) WRITE (P01REC,FMT=99992) '** MAXIT (',
     *          MAXIT, ' ) ITERATIONS HAVE BEEN PERFORMED'
            IF (IERROR.EQ.4) WRITE (P01REC,FMT=99993)
     *         '** ONE OF THE ELEMENTS OF A HAS EXCEEDED 10 IN ABSOLUTE'
     *          , ' VALUE (HEYWOOD CASE)'
            IF (IERROR.EQ.5) WRITE (P01REC,FMT=99992)
     *          '** FAILURE TO INVERT HESSIAN MATRIX AND MAXIT(', MAXIT,
     *          ') ITERATIONS MADE'
            IF (IERROR.EQ.6) WRITE (P01REC,FMT=99993)
     *          '** FAILURE TO INVERT HESSIAN MATRIX PLUS HEYWOOD CASE',
     *          ' ENCOUNTERED'
C
            IF (IERROR.EQ.7) WRITE (P01REC,FMT=99992)
     *          '** CHI-SQUARED STATISTIC HAS ', IDF,
     *          ' DEGREES OF FREEDOM'
            IFAIL = P01ABF(IEXIT,IERROR,SRNAME,1,P01REC)
C
         END IF
C
      END IF
C
      RETURN
C
99999 FORMAT (' *** ON ENTRY, THE VALUE OF CGETOL ( ',1P,D12.5,' ) IS ',
     *  'TOO SMALL',/' AND THE VALUE ',1P,D12.5,' HAS BEEN USED INS',
     *  'TEAD')
99998 FORMAT (' ** ON ENTRY, ROWS ',I16,' AND ',I16,' OF X ARE IDENTIC',
     *  'AL')
99997 FORMAT (' ** ON ENTRY, THE VALUE LW (',I16,') MUST BE AT LEAST',
     *  I16)
99996 FORMAT (' ** ON ENTRY, THE VALUE OF IRL(',I16,' ) IS ',I16)
99995 FORMAT (' ** ON ENTRY, SUM OF IRL ELEMENTS (',I16,' ) .NE. N (',
     *  I16,' )')
99994 FORMAT (' ** ON ENTRY, ONE OR MORE OF THE FOLLOWING PARAMETER VA',
     *  'LUES IS ILLEGAL',/' IP    = ',I16,'   N     = ',I16,
     *  '   IS  = ',I16,/' ISHOW = ',I16,'   NRE   = ',I16,'   ICM = ',
     *  I16,/' NRX   = ',I16,'   MAXIT = ',I16)
99993 FORMAT (1X,A,A)
99992 FORMAT (1X,A,I16,A)
99991 FORMAT (/' VALUE OF IFAIL PARAMETER ON EXIT FROM G11SAF = ',I3)
99990 FORMAT (/' ITERATION NUMBER =     0',//' THE NUMBER OF QUADRATUR',
     *  'E POINTS HAS BEEN SET TO 10')
99989 FORMAT (/' CURRENT ESTIMATES OF ALPHA(J,1)''S',/' --------------',
     *  '-------------------',/)
99988 FORMAT (' ',5F14.3)
99987 FORMAT (/' CURRENT ESTIMATES OF ALPHA(J,0)''S',/' --------------',
     *  '-------------------',/)
99986 FORMAT (/' LOG LIKELIHOOD KERNEL ON EXIT = ',D14.5)
99985 FORMAT (/' MAXIMUM LIKELIHOOD ESTIMATES OF ITEM PARAMETERS ARE A',
     *  'S FOLLOWS',/' -----------------------------------------------',
     *  '---------------',/)
99984 FORMAT (/' ITEM J',3X,'ALPHA(J,1)',5X,'ALPHA(J)',5X,'S.E.',5X,
     *  'ALPHA(J,0)',5X,'GAMMA(J)',5X,'S.E.',/' ------',3X,'----------',
     *  5X,'--------',5X,'----',5X,'----------',5X,'--------',5X,'----')
99983 FORMAT (/' ',I4,F12.3,F15.3,F10.3,F13.3,F14.3,F10.3)
99982 FORMAT (/' ITEM J',5X,'ALPHA(J)',9X,'S.E.',7X,'ALPHA(J,0)',9X,
     *  'PI(J)',10X,'S.E.',/' ------',5X,'--------',9X,'----',7X,'----',
     *  '------',9X,'-----',10X,'----')
99981 FORMAT (/' ',I4,F13.3,F15.3,F15.3,F16.3,F14.3)
      END
