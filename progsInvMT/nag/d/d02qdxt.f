      SUBROUTINE D02QDX(X,XEND,N,Y,CIN,COMM,CONST,COUT,JACSTR,JCEVAL,W,
     *                  IW,ML,MU,LRWORK,NJCPVT,NWKJAC,LWKJC0,LWKEND)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      INTEGER           MAXORD, NY2DIM, MAXSTP, MAXHNL
      CHARACTER*6       NORM, METHOD
      LOGICAL           PETZLD
      PARAMETER         (MAXORD=5,NY2DIM=6,MAXSTP=0,MAXHNL=0,
     *                  NORM='AVERAG',METHOD='NEWTON',PETZLD=.FALSE.)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, XEND
      INTEGER           IW, LRWORK, LWKEND, LWKJC0, ML, MU, N, NJCPVT,
     *                  NWKJAC
      CHARACTER*1       JACSTR
      CHARACTER*6       JCEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  CIN(7), COMM(5), CONST(6), COUT(16), W(IW), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  H0, HMAX, HMIN, TCRIT, YTEMP
      INTEGER           I, IFNEW, NEQMAX, NWLEFT
C     .. Local Arrays ..
      DOUBLE PRECISION  CONSTX(6)
C     .. External Subroutines ..
      EXTERNAL          D02NSF, D02NTF, D02NVF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN
C     .. Data statements ..
      DATA              CONSTX/2.0D0, 1.0D1, 1.0D4, 1.2D0, 1.3D0, 1.4D0/
C     .. Executable Statements ..
      DO 20 I = 1, 6
         IF (CONST(I).LT.0.0D0) THEN
            CIN(1) = -7.0D0
            GO TO 60
         ELSE IF (CONST(I).EQ.0.0D0) THEN
            CONST(I) = CONSTX(I)
         END IF
         IF (CONST(I).LE.1.0D0 .AND. I.GT.1) THEN
            CIN(1) = -8.0D0
            GO TO 60
         END IF
   20 CONTINUE
      IF (CONST(1).GE.CONST(2) .OR. CONST(2).GE.CONST(3)) THEN
         CIN(1) = -8.0D0
         GO TO 60
      END IF
      COUT(3) = 0.0D0
      YTEMP = 0.0D0
      DO 40 I = 1, N
         YTEMP = MAX(YTEMP,ABS(Y(I)))
   40 CONTINUE
      COUT(6) = YTEMP
      COUT(7) = YTEMP
      COUT(8) = 0.0D0
      COUT(9) = 0.0D0
      COUT(10) = 0.0D0
      H0 = ABS(CIN(5))*SIGN(1.0D0,XEND-X)
      TCRIT = XEND
      HMIN = COUT(13)
      HMAX = COUT(14)
      NEQMAX = N
      IFNEW = 0
      CALL D02NVF(NEQMAX,NY2DIM,MAXORD,METHOD,PETZLD,CONST,TCRIT,HMIN,
     *            HMAX,H0,MAXSTP,MAXHNL,NORM,W(LRWORK),IFNEW)
      NWLEFT = IW - LWKJC0
      IF (JACSTR.EQ.'B') THEN
         NWKJAC = (2*ML+MU+1)*N
         IF (ML.LT.0 .OR. ML.GT.N-1 .OR. MU.LT.0 .OR. MU.GT.N-1 .OR.
     *       NWKJAC.GT.NWLEFT) THEN
            CIN(1) = -1.0D0
            GO TO 60
         END IF
         IFNEW = 0
         CALL D02NTF(N,NEQMAX,JCEVAL,ML,MU,NWKJAC,NJCPVT,W(LRWORK),
     *               IFNEW)
         LWKEND = LWKJC0 + (ML+MU+1)*N
      ELSE
         NWKJAC = NEQMAX*(NEQMAX+1)
         IF (NWKJAC.GT.NWLEFT) THEN
            CIN(1) = -1.0D0
            GO TO 60
         END IF
         IFNEW = 0
         CALL D02NSF(N,NEQMAX,JCEVAL,NWKJAC,W(LRWORK),IFNEW)
         LWKEND = LWKJC0 + NEQMAX*NEQMAX
      END IF
   60 RETURN
      END
