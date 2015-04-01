      DOUBLE PRECISION FUNCTION G05EWF(R,NR,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-748 (DEC 1989).
C
C     G05EWF GENERATES THE NEXT TERM FROM AN AUTOREGRESSIVE
C     MOVING-AVERAGE TIME-SERIES USING A REFERENCE VECTOR
C     SET UP BY G05EGF.
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G05EWF')
C     .. Scalar Arguments ..
      INTEGER                          IFAIL, NR
C     .. Array Arguments ..
      DOUBLE PRECISION                 R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION                 HALF, ONE, X, ZERO
      INTEGER                          I, IPTR, J, K, MIN, MMAX, NA, NB
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G05DDF
      INTEGER                          P01ABF
      EXTERNAL                         G05DDF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        MAX, DBLE, INT
C     .. Data statements ..
      DATA                             ZERO/0.0D0/, HALF/0.5D0/,
     *                                 ONE/1.0D0/
C     .. Executable Statements ..
C
C     CHECK ARGUMENTS FOR SENSIBLE VALUES AND CONSISTENCY
C
      NA = INT(R(1))
      NB = INT(R(2))
      IPTR = INT(R(3))
      MIN = NA + NB + 5
      MMAX = MIN + MAX(NA,NB) - 1
      IF ((NA.LT.0) .OR. (NB.LT.1) .OR. (NR.LT.MMAX) .OR. (IPTR.LT.MIN)
     *     .OR. (IPTR.GT.MMAX)) GO TO 60
C
C     GENERATE THE NEXT TERM OF THE PURE AUTOREGRESSIVE SERIES.
C
      X = G05DDF(ZERO,ONE)
      J = IPTR
      DO 20 I = 1, NA
         X = X + R(I+4)*R(J)
         J = J - 1
         IF (J.LT.MIN) J = MMAX
   20 CONTINUE
      IPTR = IPTR + 1
      IF (IPTR.GT.MMAX) IPTR = MIN
      R(IPTR) = X
      R(3) = DBLE(IPTR) + HALF
C
C     NOW PERFORM A MOVING AVERAGE ON THIS.
C
      X = R(4)
      J = NA + 5
      K = J + NB - 1
      DO 40 I = J, K
         X = X + R(I)*R(IPTR)
         IPTR = IPTR - 1
         IF (IPTR.LT.MIN) IPTR = MMAX
   40 CONTINUE
      IFAIL = 0
      G05EWF = X
      RETURN
   60 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      G05EWF = ZERO
      RETURN
      END
