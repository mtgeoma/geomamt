      SUBROUTINE C06LCF(T,SIGMA,B,M,ACOEF,ERRVEC,FINV,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     C06LCF evaluates, for a specified nonnegative t, the inverse
C     finv of the prescribed Laplace transform using the series below,
C     denoted as SUM, whose coefficients are provided by C06LBF.  The
C     values of sigma and b are those returned by C06LBF.
C
C        finv = exp(sigma*t)*SUM
C
C     where
C
C        SUM = summation of (acoef(j)*exp(-b*t/2)*L(j,b*t),j=1,m),
C        L(j,b*t) denotes the Laguerre polynomial of degree j-1,
C        exp(-b*t/2)*L(j,b*t) is the associated Laguerre function.
C
C     When t is nonpositive, the evaluation approximates
C     the analytic continuation of the inverse Laplace transform,
C     becoming progressively poorer as t becomes more negative.
C
C     Note that this routine is overflow/underflow(destructive) free
C     and can be used even when the value exp(sigma*t) overflows
C     or exp(-bt/2) underflows.
C
C     C06LCF is derived from the subroutine MODUL2 in the package WEEKS
C     by B.S. Garbow, G. Giunta, J.N. Lyness and A. Murli, Algorithm
C     662: A Fortran software package for the numerical inversion of the
C     Laplace Transform based on Weeks' method, ACM Trans. Math.
C     Software, 14, pp 171-176 (1988).
C
C     INPUT arguments
C
C     t      - real - the point where the inverse Laplace transform
C                     is to be computed.
C     m      - integer - the number of terms of the Laguerre expansion.
C     acoef  - real(m) - the coefficients of the Laguerre expansion.
C     sigma  - real - the first parameter of the Laguerre expansion.
C                     It must have the same value as returned by C06LBF.
C     b      - real - the second parameter of the Laguerre expansion.
C                     It must have the same value as returned by C06LBF.
C     errvec - real(8) - the vector of diagnostic information from
C                     C06LBF. Only components 1,7 and 8 are used in
C                     C06LCF. (If C06LCF is used independently of
C                     C06LBF, store ALPHA,BETA into (7),(8) and set
C                     errvec(1) = 0.0 .)
C
C     OUTPUT arguments
C
C     finv   - real - the value of the inverse Laplace transform at t.
C     ifail  - integer - the error indicator.
C           0 => Normal termination.
C           1 => The value of the inverse Laplace transform is found to
C                be too large to be representable - finv is set to 0.0.
C           2 => The value of the inverse Laplace transform is found to
C                be too small to be representable - finv is set to 0.0.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06LCF')
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  B, FINV, SIGMA, T
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOEF(M), ERRVEC(8)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFA, BETA, BLOG, BT, BTHALF, ERRLOG, ESCALE,
     *                  EXPON, OVLOG, POLNEX, POLNOW, POLPRE, SCALE,
     *                  SLOG, SUM, THRESH, TLOG, UNLOG, UPEXP, XM
      INTEGER           IERR, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, SIGN
C     .. Executable Statements ..
C
C     Set machine-dependent constants
C
      OVLOG = -LOG(X02AMF())
      UNLOG = OVLOG
C
C     Initialize variables.
C
      ALFA = ERRVEC(7)
      BETA = ERRVEC(8)
      SCALE = UNLOG + BETA
      ESCALE = EXP(-SCALE)
      UPEXP = ZERO
      FINV = ZERO
C
      BLOG = LOG(MAX(ABS(B),ONE))
      SLOG = LOG(MAX(ABS(SIGMA),ONE))
      TLOG = LOG(MAX(ABS(T),ONE))
      IF (BLOG+TLOG.GT.OVLOG .OR. SLOG+TLOG+1.GT.OVLOG) GO TO 60
      IF (ERRVEC(1).NE.ZERO) THEN
         ERRLOG = LOG(ABS(ERRVEC(1)))
         IF (SIGMA*T+ERRLOG.GT.OVLOG) GO TO 60
         IF (SIGMA*T+ERRLOG+ALFA-BETA.LT.-UNLOG) GO TO 80
      END IF
C
C     Compute sum.
C
      BT = B*T
      BTHALF = 0.5D0*BT
      POLPRE = ZERO
      POLNOW = ONE
      SUM = ACOEF(1)
      XM = M
      THRESH = OVLOG - ALFA - LOG(MAX(XM,2+ABS(BT)))
      IF (ABS(BTHALF).LE.THRESH) THEN
         DO 20 J = 1, M - 1
            POLNEX = 2*POLNOW - POLPRE - ((1+BT)*POLNOW-POLPRE)/J
            SUM = SUM + ACOEF(J+1)*POLNEX
            POLPRE = POLNOW
            POLNOW = POLNEX
   20    CONTINUE
      ELSE
         THRESH = EXP(OVLOG-ALFA)/(2+ABS(BT))
         DO 40 J = 1, M - 1
            IF (ABS(POLNOW).GT.THRESH) THEN
               POLNOW = ESCALE*POLNOW
               POLPRE = ESCALE*POLPRE
               SUM = ESCALE*SUM
               UPEXP = UPEXP + SCALE
            END IF
            POLNEX = 2*POLNOW - POLPRE - ((1+BT)*POLNOW-POLPRE)/J
            SUM = SUM + ACOEF(J+1)*POLNEX
            POLPRE = POLNOW
            POLNOW = POLNEX
   40    CONTINUE
      END IF
C
C     Compute finv.
C
      IERR = 0
      IF (SUM.EQ.ZERO) RETURN
      EXPON = SIGMA*T - BTHALF + UPEXP + LOG(ABS(SUM))
      IF (EXPON.GT.OVLOG) GO TO 60
      IF (EXPON.LT.-UNLOG) GO TO 80
      FINV = SIGN(EXP(EXPON),SUM)
      IFAIL = 0
      GO TO 120
C
   60 IERR = 1
      WRITE (P01REC,FMT=99999)
      GO TO 100
C
   80 IERR = 2
      WRITE (P01REC,FMT=99998)
C
  100 NREC = 1
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
C
  120 RETURN
C
99999 FORMAT (' ** The approximation to f(t) is too large to be repres',
     *       'entable.')
99998 FORMAT (' ** The approximation to f(t) is too small to be repres',
     *       'entable.')
      END
