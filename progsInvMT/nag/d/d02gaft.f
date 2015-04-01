      SUBROUTINE D02GAF(U,V,N,A,B,TOL,FCN,MNP,X,Y,NP,W,LW,IW,LIW,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-305 (SEP 1981).
C     MARK 9C REVISED. IER-373 (JUN 1982)
C     MARK 10 REVISED. IER-377 (JUN 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     FCN
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02GAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, TOL
      INTEGER           IFAIL, LIW, LW, MNP, N, NP
C     .. Array Arguments ..
      DOUBLE PRECISION  U(N,2), V(N,2), W(LW), X(MNP), Y(N,MNP)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  DELEPS, H
      INTEGER           I, IND, INIT, J, LIN, LP, LWORK, M1, M2, M3, MP,
     *                  N1, NUMBEG, NUMMIX
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02GAX, D02GAY, D02GAZ, D02RAQ, D02RAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE
C     .. Executable Statements ..
      IF (N.GT.0 .AND. TOL.GT.0.0D0 .AND. (NP.GE.4 .OR. NP.EQ.0)
     *    .AND. NP.LE.MNP .AND. MNP.GE.32 .AND. LW.GE.MNP*(3*N*N+6*N+2)
     *    +4*N*N+4*N .AND. LIW.GE.MNP*(2*N+1)+N*N+4*N+2 .AND. A.LT.B)
     *    GO TO 40
   20 IND = 1
      GO TO 300
   40 NUMMIX = 0
      NUMBEG = 0
      J = 0
      DO 80 I = 1, N
         IF (V(I,1).NE.0.0D0) GO TO 60
         NUMBEG = NUMBEG + 1
   60    IF (V(I,2).NE.0.0D0) GO TO 80
         J = J + 1
   80 CONTINUE
      IF (NUMBEG.EQ.N .OR. J.EQ.N) GO TO 20
      IF (NUMBEG+J.NE.N) GO TO 20
      IF (NP.EQ.0) GO TO 120
      IF (A.NE.X(1) .OR. B.NE.X(NP)) GO TO 20
      N1 = NP - 1
      DO 100 I = 1, N1
         IF (X(I+1).LE.X(I)) GO TO 20
  100 CONTINUE
      GO TO 160
  120 NP = 4
      H = (B-A)/DBLE(NP-1)
      DO 140 I = 1, NP
         X(I) = A + DBLE(I-1)*H
  140 CONTINUE
  160 DO 200 I = 1, NP
         DO 180 J = 1, N
            Y(J,I) = ((X(I)-A)*U(J,2)+(B-X(I))*U(J,1))/(B-A)
  180    CONTINUE
  200 CONTINUE
      INIT = 1
      LIN = 2
      DELEPS = 0.0D0
      LP = MOD(IFAIL/10,10)
      MP = MOD(IFAIL/100,10)
      M1 = 1 + N
      M2 = M1 + N*N
      M3 = M2 + N*N
      LWORK = LW - 2*N*N - N
      CALL D02RAZ(N,MNP,NP,NUMBEG,NUMMIX,A,B,TOL,X,Y,N,W(1)
     *            ,FCN,D02GAZ,D02GAX,D02RAQ,D02RAQ,D02GAZ,D02GAY,D02GAZ,
     *            D02GAX,W,W,W(M1),W(M2),W,U,V,W(M3)
     *            ,LWORK,IW,LIW,DELEPS,LP,MP,INIT,LIN,IND)
      IF (IND.EQ.0) GO TO 300
      GO TO (220,240,260,280,220,220,220,220,220) IND
  220 IND = 5
      GO TO 300
  240 IND = 4
      GO TO 300
  260 IND = 2
      GO TO 300
  280 IND = 3
  300 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
