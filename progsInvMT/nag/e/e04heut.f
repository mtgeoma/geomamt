      SUBROUTINE E04HEU(IFLAG,M,N,LB,LSHESS,FVEC,X,B,IW,LIW,W,LW)
C     MARK 8 RELEASE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04HEU, AN AUXILIARY ROUTINE FOR E04GDF AND E04HEF.
C     IT EVALUATES THE SECOND DERIVATIVE TERMS IN THE HESSIAN
C     MATRIX OF THE FUNCTION.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN,
C     ENID M. R. LONG AND BRIAN T. HINDE
C     D.N.A.C.S., NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     **************************************************************
C
C     LSHESS
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LB, LIW, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LB), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSHESS
C     .. Executable Statements ..
      CALL LSHESS(IFLAG,M,N,FVEC,X,B,LB,IW,LIW,W,LW)
      RETURN
C
C     END OF E04HEU (LSHSS1)
C
      END
