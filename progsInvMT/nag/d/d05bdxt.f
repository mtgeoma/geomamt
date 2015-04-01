      SUBROUTINE D05BDX(GC,YS,N,INTS,N1,S,H,VG,VKC,WK0,FORS,TOLNL,IFLAG)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine solves the Abel equation (1). YS(i) approximates
C     the true solution  y(t), t = i*H,  i = INTS + S, INTS + S + 1, ..
C     .., INTS + S + N - 1. At each stage the value  g(i*H, YS(i)) is
C     stored in VG(i).
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, TOLNL, WK0
      INTEGER           IFLAG, INTS, N, N1, S
      CHARACTER         FORS
C     .. Array Arguments ..
      DOUBLE PRECISION  VG(0:N1+S-1), VKC(N-1), YS(0:N1+S-1)
C     .. Function Arguments ..
      DOUBLE PRECISION  GC
      EXTERNAL          GC
C     .. Local Scalars ..
      DOUBLE PRECISION  FUNYN, SCALE, SUML, TN, YSOL, YTOL
      INTEGER           IFAIL, IND, INTSS, IR, JJ, JJ1, JL
C     .. Local Arrays ..
      DOUBLE PRECISION  C(26)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          C05AXF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
      INTSS = INTS + S
      DO 60 JJ = INTSS, INTSS + N - 1
         JJ1 = JJ - 1
         TN = JJ*H
C
C        ... Start the direct summation for the remaining lag term. ...
C
         SUML = 0.0D0
C
         DO 20 JL = INTSS, JJ1
            SUML = SUML + VG(JL)*VKC(JJ-JL)
   20    CONTINUE
C
         SUML = SUML + YS(JJ)
C
C        ... Solve the nonlinear equation. ...
C
         IFAIL = 1
         SCALE = SQRT(X02AJF())
         IR = 0
         YSOL = YS(JJ-1)
         YTOL = TOLNL
         IND = 1
C
   40    CALL C05AXF(YSOL,FUNYN,YTOL,IR,SCALE,C,IND,IFAIL)
C
         IF (IND.NE.0) THEN
            IF (FORS.EQ.'S') THEN
               FUNYN = YSOL - WK0*GC(TN,YSOL) - SUML
            ELSE IF (FORS.EQ.'F') THEN
               FUNYN = -WK0*GC(TN,YSOL) - SUML
            END IF
            GO TO 40
         ELSE
            IF (IFAIL.NE.0) THEN
               IFLAG = 3
               RETURN
            END IF
         END IF
C
         VG(JJ) = GC(TN,YSOL)
         YS(JJ) = YSOL
   60 CONTINUE
C
      RETURN
      END
