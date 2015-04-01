      SUBROUTINE E04HFZ(IFLAG,M,N,FVEC,X,B,LB,IW,LIW,W,LW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     THIS SUBROUTINE COMPUTES THE SECOND TERM OF THE HESSIAN MATRIX
C     FOR THE EASY-TO-USE ROUTINE.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN
C     AND NICHOLAS I. M. GOULD
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LB, LIW, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LB), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. External Subroutines ..
      EXTERNAL          LSHES2
C     .. Executable Statements ..
      CALL LSHES2(M,N,FVEC,X,B,LB)
      RETURN
C
C     END OF E04HFZ   (TEMPHS)
C
      END
