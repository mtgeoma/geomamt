      SUBROUTINE E02DAF(M,PX,PY,X,Y,F,W,LAMDA,MU,POINT,NPOINT,DL,C,NC,
     *                  WS,NWS,EPS,SIGMA,RANK,IFAIL)
C     MARK 6 RELEASE. NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS IS A CALLING ROUTINE WHICH SETS UP WORK SPACE ARRAYS
C     FOR THE SUBROUTINE E02DAZ, AND THEN CALLS IT.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02DAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, SIGMA
      INTEGER           IFAIL, M, NC, NPOINT, NWS, PX, PY, RANK
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NC), DL(NC), F(M), LAMDA(PX), MU(PY), W(M),
     *                  WS(NWS), X(M), Y(M)
      INTEGER           POINT(NPOINT)
C     .. Local Scalars ..
      INTEGER           BW, IA, IB, ID, IDT, IL, IRT, IT, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02DAZ
C     .. Executable Statements ..
      BW = 4 + 3*(PY-4)
      K = 0
      IF (NWS.LT.2*NC*(2+BW)+BW) K = 4
      IF (K.GT.0) GO TO 20
      ID = 1
      IB = ID + NC
      IT = IB + NC
      IA = IT + NC*BW
      IDT = IA + BW
      IL = IDT + NC
      IRT = IL + NC
      CALL E02DAZ(M,PX,PY,X,Y,F,W,LAMDA,MU,POINT,NPOINT,DL,C,NC,WS(ID)
     *            ,WS(IB),WS(IT),WS(IA),BW,WS(IDT),WS(IL),WS(IRT)
     *            ,EPS,SIGMA,RANK,K)
      IF (K.NE.0) GO TO 20
      IFAIL = K
      RETURN
   20 IFAIL = P01ABF(IFAIL,K,SRNAME,0,P01REC)
      RETURN
      END
