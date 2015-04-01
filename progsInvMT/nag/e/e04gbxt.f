      SUBROUTINE E04GBX(IFLAG,M,N,LSFJAC,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ************************************************************
C
C     E04GBX, AN AUXILIARY ROUTINE FOR SUBROUTINES E04GBF, E04GDF,
C     AND E04HEF EVALUATES THE FUNCTION VECTOR AND JACOBIAN
C     MATRIX.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN
C     ENID M. R. LONG AND BRIAN T. HINDE
C     D.N.A.C.S., NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     ************************************************************
C
C     LSFJAC
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LIW, LJ, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSFJAC
C     .. Executable Statements ..
      CALL LSFJAC(IFLAG,M,N,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
      RETURN
C
C     END OF E04GBX  (LSFJS1)
C
      END
