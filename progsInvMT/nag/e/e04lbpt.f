      SUBROUTINE E04LBP(IFLAG,N,LH,X,G,ISTATE,SFUN,DELTA,HESL,HESD,IW,
     *                  LIW,W,LW)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBP (ZHESS) IS CALLED FROM E04LBR TO COMPUTE A
C     FINITE-DIFFERENCE APPROXIMATION OF THE PROJECTED HESSIAN
C     MATRIX AND OF SUCH ADDITIONAL ELEMENTS OF THE FULL HESSIAN AS
C     WILL BE NEEDED FOR THE CALCULATION OF SECOND-ORDER LAGRANGE
C     MULTIPLIERS.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     SFUN
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DELTA
      INTEGER           IFLAG, LH, LIW, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), HESD(N), HESL(LH), W(LW), X(N)
      INTEGER           ISTATE(N), IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          SFUN
C     .. Local Scalars ..
      DOUBLE PRECISION  DINV, DINV2, GI, XI
      INTEGER           I, IH, IMINUS, IPLUS, ISJ, ISTART, J, JUMP
C     .. Executable Statements ..
      DO 20 I = 1, N
         HESD(I) = 0.0D+0
   20 CONTINUE
      DO 40 I = 1, LH
         HESL(I) = 0.0D+0
   40 CONTINUE
      DINV = 1.0D+0/DELTA
      DINV2 = 0.5D+0*DINV
      ISTART = 0
      DO 180 I = 1, N
         ISTART = ISTART + I
         IF (ISTATE(I).LT.0) GO TO 180
         XI = X(I)
         X(I) = XI + DELTA
         IFLAG = 1
         IW(2) = I
         CALL SFUN(IFLAG,N,X,GI,W,IW,LIW,W,LW)
         IF (IFLAG.LT.0) RETURN
         X(I) = XI
         GI = G(I)
         HESD(I) = (W(I)-GI)*DINV
         IF (I.EQ.N) GO TO 120
         IH = ISTART
         JUMP = I
         IPLUS = I + 1
         DO 100 J = IPLUS, N
            ISJ = ISTATE(J)
            IF (ISJ.EQ.-3) GO TO 80
            IF (ISJ.LT.0) GO TO 60
            HESL(IH) = W(J) - GI
            GO TO 80
   60       HESL(IH) = (W(J)-G(J))*DINV
   80       IH = IH + JUMP
            JUMP = JUMP + 1
  100    CONTINUE
  120    IF (I.EQ.1) GO TO 180
         IMINUS = I - 1
         IH = IMINUS*(IMINUS-1)/2
         DO 160 J = 1, IMINUS
            IH = IH + 1
            ISJ = ISTATE(J)
            IF (ISJ.EQ.-3) GO TO 160
            IF (ISJ.LT.0) GO TO 140
            HESL(IH) = (HESL(IH)+W(J)-GI)*DINV2
            GO TO 160
  140       HESL(IH) = (W(J)-G(J))*DINV
  160    CONTINUE
  180 CONTINUE
      RETURN
C
C     END OF E04LBP (ZHESS)
C
      END
