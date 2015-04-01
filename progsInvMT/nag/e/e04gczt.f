      SUBROUTINE E04GCZ(IFLAG,M,N,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     THIS SUBROUTINE COMPUTES THE FUNCTION VECTOR AND JACOBIAN
C     MATRIX FOR THE EASY-TO-USE ROUTINES.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN
C     AND NICHOLAS I. M. GOULD
C     D.N.A.C.,NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LIW, LJ, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. External Subroutines ..
      EXTERNAL          LSFUN2
C     .. Executable Statements ..
      CALL LSFUN2(M,N,X,FVEC,FJAC,LJ)
      RETURN
C
C     END OF E04GCZ   (TEMPFN)
C
      END
