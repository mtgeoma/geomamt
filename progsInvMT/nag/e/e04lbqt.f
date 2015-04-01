      SUBROUTINE E04LBQ(IFLAG,N,LH,X,G,ISTATE,INEW,SFUN,DELTA,HESL,HESD,
     *                  IW,LIW,W,LW)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBQ (UPZHSS) IS CALLED FROM E04LBR TO UPDATE THE
C     FINITE-DIFFERENCE APPROXIMATION OF THE PROJECTED HESSIAN
C     MATRIX, AND TO COMPUTE SUCH ADDITIONAL ELEMENTS OF THE FULL
C     HESSIAN AS WILL BE NEEDED FOR THE CALCULATION OF SECOND-ORDER
C     LAGRANGE MULTIPLIERS, WHEN A FIXED VARIABLE X(INEW) HAS BEEN
C     RELEASED FROM ITS BOUND BUT THE VALUES OF THE FREE VARIABLES
C     HAVE NOT CHANGED.
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
      INTEGER           IFLAG, INEW, LH, LIW, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), HESD(N), HESL(LH), W(LW), X(N)
      INTEGER           ISTATE(N), IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          SFUN
C     .. Local Scalars ..
      DOUBLE PRECISION  DINV, TEMP, XNEW
      INTEGER           IH, IMINUS, IPLUS, ISJ, J, JUMP
C     .. Executable Statements ..
      DINV = 1.0D+0/DELTA
      XNEW = X(INEW)
      X(INEW) = XNEW + DELTA
      IFLAG = 1
      IW(2) = INEW
      CALL SFUN(IFLAG,N,X,TEMP,W,IW,LIW,W,LW)
      IF (IFLAG.LT.0) RETURN
      X(INEW) = XNEW
      HESD(INEW) = (W(INEW)-G(INEW))*DINV
      IF (INEW.EQ.1) GO TO 60
      IMINUS = INEW - 1
      IH = IMINUS*(IMINUS-1)/2
      DO 40 J = 1, IMINUS
         IH = IH + 1
         ISJ = ISTATE(J)
         IF (ISJ.EQ.-3) GO TO 40
         TEMP = (W(J)-G(J))*DINV
         IF (ISJ.LT.0) GO TO 20
         HESL(IH) = 0.5D+0*(HESL(IH)+TEMP)
         GO TO 40
   20    HESL(IH) = TEMP
   40 CONTINUE
   60 IF (INEW.EQ.N) RETURN
      IPLUS = INEW + 1
      IH = INEW*IPLUS/2
      JUMP = INEW
      DO 120 J = IPLUS, N
         ISJ = ISTATE(J)
         IF (ISJ.EQ.-3) GO TO 100
         TEMP = (W(J)-G(J))*DINV
         IF (ISJ.LT.0) GO TO 80
         HESL(IH) = 0.5D+0*(HESL(IH)+TEMP)
         GO TO 100
   80    HESL(IH) = TEMP
  100    IH = IH + JUMP
         JUMP = JUMP + 1
  120 CONTINUE
      RETURN
C
C     END OF E04LBQ (UPZHSS)
C
      END
