      SUBROUTINE D02RAF(M,NMAX,N,NUMBEG,NUMMIX,TOL,INIT,X,Y,IY,ABT,FCN,
     *                  G,IJAC,JACOBF,JACOBG,DELEPS,JACEPS,JACGEP,WORK,
     *                  LWORK,IWORK,LIWORK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-821 (DEC 1989).
C     FCN, G, JACEPS, JACGEP, JACOBF, JACOBG
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02RAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DELEPS, TOL
      INTEGER           IFAIL, IJAC, INIT, IY, LIWORK, LWORK, M, N,
     *                  NMAX, NUMBEG, NUMMIX
C     .. Array Arguments ..
      DOUBLE PRECISION  ABT(M), WORK(LWORK), X(NMAX), Y(IY,NMAX)
      INTEGER           IWORK(LIWORK)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, G, JACEPS, JACGEP, JACOBF, JACOBG
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B
      INTEGER           I, IFLAG, LIN, LP, LW, M1, M2, MP, N1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02PAZ, X02AJF
      INTEGER           P01ABF
      EXTERNAL          D02PAZ, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02HBS, D02RAQ, D02RAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, SQRT
C     .. Executable Statements ..
      IF ((TOL.LE.0.0D0) .OR. (IY.LT.M)) GO TO 60
      A = D02PAZ(SQRT(X02AJF()))
      IF (DELEPS.GT.0.0D0 .AND. DELEPS.LT.A)
     *    GO TO 60
      IF (DELEPS.LT.0.0D0 .OR. DELEPS.GE.1.0D0) DELEPS = 0.0D0
      LIN = 4
      IF (IJAC.NE.0) LIN = 3
      LP = MOD(IFAIL/10,10)
      MP = MOD(IFAIL/100,10)
      IF (X(1).GE.X(N)) GO TO 60
      IF (INIT.EQ.0) GO TO 40
      N1 = N - 1
      DO 20 I = 1, N1
         IF (X(I+1).LE.X(I)) GO TO 60
   20 CONTINUE
   40 A = X(1)
      B = X(N)
      LW = LWORK - 2*M*M
      M1 = 1 + M*M
      M2 = M1 + M*M
      CALL D02RAZ(M,NMAX,N,NUMBEG,NUMMIX,A,B,TOL,X,Y,IY,ABT,D02HBS,G,
     *            FCN,D02RAQ,D02RAQ,JACOBF,JACOBG,JACEPS,JACGEP,WORK,
     *            WORK,WORK(1),WORK(M1),WORK,WORK,WORK,WORK(M2)
     *            ,LW,IWORK,LIWORK,DELEPS,LP,MP,INIT,LIN,IFLAG)
      IF (IFLAG.EQ.0) GO TO 160
      GO TO (60,160,160,160,60,80,100,120,140)
     *       IFLAG
   60 IFLAG = 1
      GO TO 160
   80 IFLAG = 9
      GO TO 160
  100 IFLAG = 6
      IF (DELEPS.EQ.0.0D0) IFLAG = 8
      GO TO 160
  120 IFLAG = 5
      GO TO 160
  140 IFLAG = 7
      IF (DELEPS.EQ.0.0D0) IFLAG = 8
  160 IFAIL = P01ABF(IFAIL,IFLAG,SRNAME,0,P01REC)
      RETURN
      END
