      SUBROUTINE E04HDF(N,SFUN,SHESS,X,G,HESL,LH,HESD,IW,LIW,W,LW,IFAIL)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 13 REVISED. IER-645 (APR 1988).
C     MARK 14 REVISED. IER-800 (DEC 1989).
C
C     **************************************************************
C
C     E04HDF CHECKS THAT A USER-SUPPLIED ROUTINE FOR CALCULATING
C     SECOND DERIVATIVES OF AN OBJECTIVE FUNCTION IS CONSISTENT WITH
C     A USER-SUPPLIED ROUTINE FOR CALCULATING THE CORRESPONDING
C     FIRST DERIVATIVES.
C
C     THE ROUTINE IS ESSENTIALLY IDENTICAL TO THE SUBROUTINE CHKHES
C     IN THE NPL ALGORITHMS LIBRARY (REF. NO. E4/04/F). W(I), I = 1,
C     2, . . . , 5*N ARE USED AS WORKSPACE. (NOTE THAT, FOR
C     CONSISTENCY WITH OTHER E04 DOCUMENTATION, THE NAMES FUNCT AND
C     HESS ARE USED INSTEAD OF SFUN AND SHESS IN THE WRITE-UP.)
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY, SUSAN M.
C     PICKEN, D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     SFUN, SHESS
C
C     A MACHINE DEPENDENT CONSTANT IS SET HERE. EPSMCH IS THE
C     SMALLEST POSITIVE REAL NUMBER SUCH THAT 1.0 + EPSMCH .GT. 1.0.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04HDF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LH, LIW, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), HESD(N), HESL(LH), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          SFUN, SHESS
C     .. Local Scalars ..
      DOUBLE PRECISION  APPRX, DEL, EPSMCH, FDUMMY, RTDEL, TEMP, YTG,
     *                  YTGY, YTHY, ZTG
      INTEGER           I, IG, IP, IWK, IXYI, IXYP, IY, IYI, IYP, IZ, K,
     *                  NWHY
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DDOT, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04HCZ, E04HDZ, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      EPSMCH = X02AJF()
      DEL = SQRT(EPSMCH)
      RTDEL = SQRT(DEL)
      IF (N.GE.1 .AND. LH.GE.N*(N-1)/2 .AND. LH.GT.0 .AND. LW.GE.5*
     *    N .AND. LIW.GT.0) GO TO 20
      NWHY = 1
      GO TO 100
C
C     THE HESSIAN MATRIX, WHICH IS TO BE CHECKED, AND THE GRADIENT
C     VECTOR G ARE CALCULATED AT X.
C
   20 NWHY = 2
      CALL SFUN(NWHY,N,X,FDUMMY,G,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 100
      CALL DCOPY(N,G,1,HESD,1)
      CALL SHESS(NWHY,N,X,HESL,LH,HESD,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 100
      IY = 1
      IZ = IY + N
      IXYP = IZ + N
      IWK = IXYP + N
      IG = IWK + N
C
C     TWO ORTHOGONAL UNIT VECTORS Y AND Z ARE NOW SET UP AND THEIR
C     INNER PRODUCTS WITH THE GRADIENT VECTOR ARE FORMED.
C
      CALL E04HCZ(N,W(IY),W(IZ))
      YTG = DDOT(N,W(IY),1,G,1)
      ZTG = DDOT(N,W(IZ),1,G,1)
C
C     PERFORM THE LOOP FOR Y THEN Z.
C
      IP = 0
      DO 80 K = 1, 2
         IF (N.EQ.1) GO TO 60
C
C        FORM XY = X + DEL*Y.
C
         IYP = IY + IP
         DO 40 I = 1, N
            IXYI = IXYP + I - 1
            IYI = IYP + I - 1
            W(IXYI) = X(I) + DEL*W(IYI)
   40    CONTINUE
C
C        FORM YT*H*Y, WHERE H IS THE HESSIAN MATRIX.
C
         CALL E04HDZ(N,LH,HESL,HESD,W(IYP),W(IWK),YTHY)
C
C        EVALUATE THE GRADIENT AT XY AND FORM ITS INNER PRODUCT WITH Y.
C
         NWHY = 2
         CALL SFUN(NWHY,N,W(IXYP),TEMP,W(IG),IW,LIW,W,LW)
         IF (NWHY.LT.0) GO TO 100
         YTGY = DDOT(N,W(IYP),1,W(IG),1)
         NWHY = 2
C
C        COMPARE THE HESSIAN MATRIX WITH A FINITE DIFFERENCE
C        APPROXIMATION AND SET NWHY = 2 IF THE HESSIAN IS INCORRECT.
C
         APPRX = (YTGY-YTG)/DEL
         IF (ABS(YTHY-APPRX).GT.RTDEL*(ABS(YTHY)+1.0D+0)) GO TO 100
   60    IP = IP + N
         YTG = ZTG
   80 CONTINUE
      NWHY = 0
  100 IF (NWHY.NE.0) GO TO 120
      IFAIL = 0
      RETURN
  120 IFAIL = P01ABF(IFAIL,NWHY,SRNAME,0,P01REC)
      RETURN
C
C     END OF E04HDF (CHKHES)
C
      END
