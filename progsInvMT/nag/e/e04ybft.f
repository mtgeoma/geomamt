      SUBROUTINE E04YBF(M,N,LSQFUN,LSQHES,X,FVEC,FJAC,LJ,B,LB,IW,LIW,W,
     *                  LW,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED.  IER-683 (DEC 1989).
C
C     **************************************************************
C
C     E04YBF CHECKS THAT A USER-SUPPLIED ROUTINE FOR EVALUATING
C     THE SECOND DERIVATIVE TERM OF THE HESSIAN MATRIX OF A SUM OF
C     SQUARES IS CONSISTENT WITH A USER-SUPPLIED ROUTINE FOR
C     CALCULATING THE CORRESPONDING FIRST DERIVATIVES.
C
C     THIS SUBROUTINE MAY BE USED WHEN THE PARAMETER LIST OF LSQHES
C     IS  (IFLAG, M, N, LB, FVEC, X, B, IW, LIW, W, LW)
C     AND THE PARAMETER LIST OF LSQFUN IS
C     (IFLAG, M, N, X, FVEC, FJAC, LJ, IW, LIW, W, LW)
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
      PARAMETER         (SRNAME='E04YBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LB, LIW, LJ, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LB), FJAC(LJ,N), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSQFUN, LSQHES
C     .. Local Scalars ..
      DOUBLE PRECISION  APPRX, DEL, EPSMCH, SUM, YTG, YTGY, YTHY, ZTG
      INTEGER           I, IFJAC, IFVEC, IG, IHESD, IHESL, IP, IXYI,
     *                  IXYP, IXYP1, IY, IYI, IYP, IYP1, IZ, J, L, LB1,
     *                  LH, LI, LL, LW1, LW2, LWORK, NREC, NWHY
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DDOT, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04HCZ, E04HDZ, E04HEZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      EPSMCH = X02AJF()
      DEL = SQRT(EPSMCH)
      NWHY = 1
      LWORK = 5*N + N*(N-1)/2 + M + M*N
      IF (N.EQ.1) LWORK = LWORK + 1
      IF (M.LT.N .OR. N.LT.1 .OR. LJ.LT.M .OR. LB.LT.N*(N+1)
     *    /2 .OR. LW.LT.LWORK .OR. LIW.LE.0) GO TO 140
C
C     SET UP THE WORK SPACE.
C
      IY = 1
      IZ = IY + N
      IXYP = IZ + N
      IG = IXYP + N
      IFVEC = IG + N
      IFJAC = IFVEC + M
      IHESD = IFJAC + M*N
      IHESL = IHESD + N
      LH = N*(N-1)/2
      IF (N.EQ.1) LH = 1
C
C     COMPUTE THE GRADIENT OF THE SUM OF SQUARES.
C
      NWHY = 2
      CALL LSQFUN(NWHY,M,N,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 140
      CALL E04HEZ(M,N,FVEC,FJAC,LJ,W(IG))
      NWHY = 0
      CALL LSQHES(NWHY,M,N,FVEC,X,B,LB,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 140
C
C     COMPUTE THE HESSIAN MATRIX OF THE SUM OF SQUARES, 2(JTJ + B),
C     AND STORE IN HESL AND HESD.
C
      L = IHESL
      LI = IHESD
      LL = 1
      DO 60 I = 1, N
         DO 40 J = 1, I
            SUM = 2.0D+0*(DDOT(M,FJAC(1,I),1,FJAC(1,J),1)+B(LL))
            LL = LL + 1
            IF (I.EQ.J) GO TO 20
            W(L) = SUM
            L = L + 1
            GO TO 40
   20       W(LI) = SUM
            LI = LI + 1
   40    CONTINUE
   60 CONTINUE
C
C     TWO ORTHOGONAL UNIT VECTORS Y AND Z ARE NOW SET UP AND THEIR
C     INNER PRODUCTS WITH THE GRADIENT VECTOR ARE FORMED.
C
      CALL E04HCZ(N,W(IY),W(IZ))
      YTG = DDOT(N,W(IY),1,W(IG),1)
      ZTG = DDOT(N,W(IZ),1,W(IG),1)
C
C     PERFORM THE LOOP FOR Y THEN Z.
C
      IP = 0
      DO 120 L = 1, 2
         IF (N.EQ.1) GO TO 100
C
C        FORM XY = X + DEL*Y.
C
         IYP = IY + IP
         IXYP1 = IXYP - 1
         IYP1 = IYP - 1
         DO 80 I = 1, N
            IXYI = IXYP1 + I
            IYI = IYP1 + I
            W(IXYI) = X(I) + DEL*W(IYI)
   80    CONTINUE
C
C        FORM YT*H*Y, WHERE H IS THE HESSIAN MATRIX.
C
         CALL E04HDZ(N,LH,W(IHESL),W(IHESD),W(IYP),W(IG),YTHY)
C
C        EVALUATE THE GRADIENT AT XY AND FORM ITS INNER PRODUCT WITH Y.
C
         NWHY = 2
         CALL LSQFUN(NWHY,M,N,W(IXYP),W(IFVEC),W(IFJAC),M,IW,LIW,W,LW)
         IF (NWHY.LT.0) GO TO 140
         CALL E04HEZ(M,N,W(IFVEC),W(IFJAC),M,W(IG))
         YTGY = DDOT(N,W(IYP),1,W(IG),1)
         NWHY = 2
C
C        COMPARE THE HESSIAN MATRIX WITH A FINITE DIFFERENCE
C        APPROXIMATION AND SET NWHY = 2 IF THE HESSIAN IS INCORRECT.
C
         APPRX = (YTGY-YTG)/DEL
         IF (ABS(YTHY-APPRX).GT.SQRT(DEL)*(ABS(YTHY)+1.0D+0))
     *       GO TO 140
  100    IP = IP + N
         YTG = ZTG
  120 CONTINUE
      NWHY = 0
  140 CONTINUE
      IF (NWHY.LT.0) THEN
         P01REC(1) =
     *     ' ** Negative value of IFLAG set in LSQFUN or LSQHES by user'
         NREC = 1
      ELSE IF (NWHY.EQ.1) THEN
         LW1 = 5*N + M + M*N + N*(N-1)/2
         LW2 = 6 + 2*M
         LB1 = (N+1)*N/2
         IF (M.LT.N) THEN
            WRITE (P01REC(1),FMT=99997) M, N
            NREC = 1
         ELSE IF (N.LT.1) THEN
            WRITE (P01REC(1),FMT=99996) N
            NREC = 1
         ELSE IF (LJ.LT.M) THEN
            WRITE (P01REC(1),FMT=99995) LJ, M
            NREC = 1
         ELSE IF (LB.LT.LB1) THEN
            WRITE (P01REC(1),FMT=99993) LB, LB1
            NREC = 1
         ELSE IF (LIW.LT.1) THEN
            WRITE (P01REC(1),FMT=99994) LIW
            NREC = 1
         ELSE IF (LW.LT.LW1 .AND. N.GT.1) THEN
            WRITE (P01REC,FMT=99999) LW, LW1
            NREC = 2
         ELSE IF (LW.LT.LW2 .AND. N.EQ.1) THEN
            WRITE (P01REC,FMT=99998) LW, LW2
            NREC = 2
         END IF
      ELSE IF (NWHY.EQ.2) THEN
         P01REC(1) =
     *        ' ** It is very likely that the user has made an error in'
         P01REC(2) = ' ** setting up the array B in LSQHES'
         NREC = 2
      END IF
      IFAIL = P01ABF(IFAIL,NWHY,SRNAME,NREC,P01REC)
      RETURN
C
C     END OF E04YBF   (CHKLSH)
C
99999 FORMAT (' ** On entry, LW must be at least 5*N + M + M*N + N*(N-',
     *       '1)/2 if N.gt.1:',/' ** LW =',I16,', LW must be at least',
     *       I16)
99998 FORMAT (' ** On entry, LW must be at least 6 + 2*M:',/' ** LW =',
     *       I16,', LW must be at least',I16)
99997 FORMAT (' ** On entry, M must be at least N: M =',I16,', N =',I16)
99996 FORMAT (' ** On entry, N must be at least 1: N =',I16)
99995 FORMAT (' ** On entry, LJ must be at least M: LJ =',I16,', M =',
     *       I16)
99994 FORMAT (' ** On entry, LIW must be at least 1: LIW =',I16,
     *       ', N = ',I16)
99993 FORMAT (' ** On entry, LB must be at least (N+1)*N/2: LB =',I16,
     *       ', (N+1)*N/2 =',I16)
      END
