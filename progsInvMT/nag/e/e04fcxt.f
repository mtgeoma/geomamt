      SUBROUTINE E04FCX(IFLAG,M,N,GETFVC,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 9 REVISED. IER-318 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     **************************************************************
C
C     E04FCX FORMS A FINITE-DIFFERENCE APPROXIMATION TO THE JACOBIAN
C     MATRIX OF FIRST DERIVATIVES OF THE FUNCTION VECTOR.
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY
C     AND SUSAN M. PICKEN
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     **************************************************************
C
C     GETFVC
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LIW, LJ, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          GETFVC
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSMCH, H, HINV, ONE, RTEPS, XJ
      INTEGER           I, J
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      ONE = 1.0D+0
      EPSMCH = X02AJF()
      RTEPS = SQRT(EPSMCH)
      DO 40 J = 1, N
         XJ = X(J)
         H = RTEPS*(ONE+ABS(XJ))
         HINV = ONE/H
         X(J) = XJ + H
         CALL GETFVC(IFLAG,M,N,X,W,IW,LIW,W,LW)
         IF (IFLAG.LT.0) RETURN
         DO 20 I = 1, M
            FJAC(I,J) = (W(I)-FVEC(I))*HINV
   20    CONTINUE
         X(J) = XJ
   40 CONTINUE
      RETURN
C
C     END OF E04FCX   (APJCBN)
C
      END
