      SUBROUTINE E04YAF(M,N,LSQFUN,X,FVEC,FJAC,LJ,IW,LIW,W,LW,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     **************************************************************
C
C     E04YAF CHECKS THAT A USER-SUPPLIED ROUTINE FOR EVALUATING A
C     VECTOR OF FUNCTIONS AND THE MATRIX OF THEIR FIRST DERIVATIVES
C     PRODUCES DERIVATIVE VALUES WHICH ARE CONSISTENT WITH THE
C     FUNCTION VALUES CALCULATED.
C
C     THIS SUBROUTINE MAY BE USED WHEN THE PARAMETER LIST OF LSQFUN
C     IS  (IFLAG, M, N, X, FVEC, FJAC, LJ, IW, LIW, W, LW)
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY,
C     SUSAN M. PICKEN
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     Revised to output explanatory messages.
C     Peter Mayes, NAG Central Office, December 1987.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04YAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LIW, LJ, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSQFUN
C     .. Local Scalars ..
      DOUBLE PRECISION  C1, C2, EPSMCH, F, F1, H, V1, V2, W1, W2
      INTEGER           I, IA, IA1, IAI, IB, IB1, IBI, IFJAC, IFVEC, IG,
     *                  LW1, NREC, NWHY
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DDOT, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04HCZ, E04HEZ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      EPSMCH = X02AJF()
      H = SQRT(EPSMCH)
      NWHY = 1
      IF (N.LT.1 .OR. M.LT.N .OR. LJ.LT.M .OR. LIW.LE.0 .OR. LW.LT.3*N+
     *    M+M*N) GO TO 80
C
C     THE SUM OF SQUARES AND GRADIENT VECTOR, EVALUATED AT X,
C     ARE ASSIGNED TO F AND ARRAY W(IG) RESPECTIVELY
C
      IA = 1
      IB = IA + N
      IG = IB + N
      IFVEC = IG + N
      IFJAC = IFVEC + M
      NWHY = 2
      CALL LSQFUN(NWHY,M,N,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 80
      CALL E04HEZ(M,N,FVEC,FJAC,LJ,W(IG))
      F = DDOT(M,FVEC,1,FVEC,1)
C
C     TWO ORTHOGONAL VECTORS A AND B ARE SET UP.
C
      CALL E04HCZ(N,W(IA),W(IB))
C
C     W1 AND W2 ARE PUT EQUAL TO GT*A AND GT*B RESPECTIVELY,
C     WHERE T DENOTES TRANSPOSE
C
      W1 = DDOT(N,W(IA),1,W(IG),1)
      W2 = DDOT(N,W(IB),1,W(IG),1)
C
C     A FORWARD-DIFFERENCE APPROXIMATION TO THE GRADIENT IS
C     NOW MADE ALONG A AND B.
C
      IA1 = IA - 1
      DO 20 I = 1, N
         IAI = IA1 + I
         W(IAI) = X(I) + H*W(IAI)
   20 CONTINUE
      F1 = F
      IF (N.EQ.1) GO TO 40
      NWHY = 2
      CALL LSQFUN(NWHY,M,N,W(IA),W(IFVEC),W(IFJAC),M,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 80
      F1 = DDOT(M,W(IFVEC),1,W(IFVEC),1)
   40 V1 = (F1-F)/H
      IB1 = IB - 1
      DO 60 I = 1, N
         IBI = IB1 + I
         W(IBI) = X(I) + H*W(IBI)
   60 CONTINUE
      NWHY = 2
      CALL LSQFUN(NWHY,M,N,W(IB),W(IFVEC),W(IFJAC),M,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 80
      F1 = DDOT(M,W(IFVEC),1,W(IFVEC),1)
      V2 = (F1-F)/H
C
C     C1 AND C2 ARE THE DIFFERENCES BETWEEN APPROXIMATED AND
C     PROGRAMMED GRADIENT PROJECTED ALONG A AND B RESPECTIVELY
C
      C1 = V1 - W1
      C2 = V2 - W2
C
C     IF EITHER C1 OR C2 IS TOO LARGE AN ERROR INDICATOR IS SET
C
      NWHY = 0
      IF (C1*C1.GE.H*(W1*W1+1.0D+0) .OR. C2*C2.GE.H*(W2*W2+1.0D+0))
     *    NWHY = 2
   80 CONTINUE
      IF (NWHY.LT.0) THEN
         P01REC(1) = ' ** Negative value of IFLAG set in LSQFUN by user'
         NREC = 1
      ELSE IF (NWHY.EQ.1) THEN
         LW1 = 3*N + M + M*N
         IF (M.LT.N) THEN
            WRITE (P01REC(1),FMT=99999) M, N
            NREC = 1
         ELSE IF (N.LT.1) THEN
            WRITE (P01REC(1),FMT=99998) N
            NREC = 1
         ELSE IF (LJ.LT.M) THEN
            WRITE (P01REC(1),FMT=99997) LJ, M
            NREC = 1
         ELSE IF (LIW.LT.1) THEN
            WRITE (P01REC(1),FMT=99996) LIW
            NREC = 1
         ELSE IF (LW.LT.LW1) THEN
            WRITE (P01REC,FMT=99995) LW, LW1
            NREC = 2
         END IF
      ELSE IF (NWHY.EQ.2) THEN
         P01REC(1) =
     *        ' ** It is very likely that the user has made an error in'
         P01REC(2) = ' ** forming the derivatives in LSQFUN'
         NREC = 2
      END IF
      IFAIL = P01ABF(IFAIL,NWHY,SRNAME,NREC,P01REC)
      RETURN
C
C     END OF E04YAF   (CHKLSJ)
C
99999 FORMAT (' ** On entry, M must be at least N: M =',I16,', N =',I16)
99998 FORMAT (' ** On entry, N must be at least 1: N =',I16)
99997 FORMAT (' ** On entry, LJ must be at least M: LJ =',I16,', M =',
     *       I16)
99996 FORMAT (' ** On entry, LIW must be at least 1: LIW =',I16,
     *       ', N = ',I16)
99995 FORMAT (' ** On entry, LW must be at least 3*N + M + M*N:',/' **',
     *       ' LW =',I16,', LW must be at least',I16)
      END
