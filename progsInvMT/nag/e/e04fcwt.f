      SUBROUTINE E04FCW(IFLAG,M,N,LSFUN,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-233 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04FCW, AN AUXILIARY ROUTINE FOR SUBROUTINE E04FCF AND E04FDF
C     EVALUATES THE FUNCTION VECTOR AND, IF NECESSARY, FORMS A
C     FINITE-DIFFERENCE APPROXIMATION TO ITS JACOBIAN MATRIX.
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY,
C     SUSAN M. PICKEN AND BRIAN T. HINDE
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     **************************************************************
C
C
C     LSFUN
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LIW, LJ, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSFUN
C     .. External Subroutines ..
      EXTERNAL          E04FCX
C     .. Executable Statements ..
      IF (IFLAG.NE.1) CALL LSFUN(IFLAG,M,N,X,FVEC,IW,LIW,W,LW)
      IF (IFLAG.LE.0) RETURN
      CALL E04FCX(IFLAG,M,N,LSFUN,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
      RETURN
C
C     END OF E04FCW   (LSFJS2)
C
      END
