      SUBROUTINE E04DGY(N,HNEW,HOLD,ALPHA,PK,YK,G,ITER,DEBUG,IDBGCG)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16A REVISED. IER-986 (JUN 1993).
C
C     E04DGY updates the scaling.
C
C     The ideas for this routine are discussed on pages 18-19
C     of the CG paper, in particular formula ( 4.4 ). omega
C     corresponds to BOUNDK. delta-bar and delta correspond
C     to HNEW and HOLD, which in turn are actually both
C     equal to DIAGB when this routine is called.
C
C
C     -- Written on 4-June-1986.
C     Sven Hammarling and Janet Welding.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA
      INTEGER           IDBGCG, ITER, N
      LOGICAL           DEBUG
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), HNEW(N), HOLD(N), PK(N), YK(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, BOUNDK, COND, EPSMCH, GI, GTP, HII, PWR,
     *                  SMALL, TINY, YKI, YTP
      INTEGER           I, J
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, LOG, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
      BOUNDK = 1.0D-2/(SQRT(DBLE(N))*EPSMCH)
      SMALL = BOUNDK
      BIG = ZERO
      YTP = DDOT(N,YK,1,PK,1)
      GTP = DDOT(N,G,1,PK,1)
      DO 20 I = 1, N
         YKI = YK(I)
         GI = G(I)
         HII = HOLD(I) + (GI**2)/GTP + (YKI**2)/(ALPHA*YTP)
         IF (HII.LT.ZERO) HII = HOLD(I)
         IF (HII.GT.BIG) BIG = HII
         IF (HII.LT.SMALL) SMALL = HII
         HNEW(I) = HII
   20 CONTINUE
      TINY = EPSMCH*BIG
      IF (SMALL.LE.TINY) THEN
C
C        Zero diagonal element ( caused by rank-two update
C        giving a negative element ).
C
         SMALL = TINY
         DO 40 I = 1, N
            IF (HNEW(I).LT.TINY) HNEW(I) = TINY
   40    CONTINUE
      END IF
      COND = BIG/SMALL
      IF (COND.LE.BOUNDK) GO TO 80
C
C     Raise each diagonal element to the power PWR where PWR is such
C     that the new ratio of the largest element to the smallest of
C     HNEW is BOUNDK.
C
      PWR = LOG(BOUNDK)/LOG(COND)
      DO 60 I = 1, N
         HNEW(I) = HNEW(I)**PWR
   60 CONTINUE
   80 CONTINUE
      IF (DEBUG .AND. ITER.GE.IDBGCG) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(NOUT,2,REC)
         DO 100 I = 1, N, 5
            WRITE (REC,FMT=99998) (HNEW(J),J=1,MIN(I+4,N))
            CALL X04BAF(NOUT,REC(1))
  100    CONTINUE
      END IF
      RETURN
C
C     End of E04DGY. (UPDATE).
C
99999 FORMAT (/' //E04DGY// - N element vector HNEW is')
99998 FORMAT (1X,5G15.7)
      END
