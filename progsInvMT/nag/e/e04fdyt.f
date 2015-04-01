      SUBROUTINE E04FDY(IFLAG,M,N,X,FVEC,IW,LIW,W,LW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     CALLS THE USER-SUPPLIED ROUTINE LSFUN1.
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY AND
C     SUSAN M. PICKEN
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LIW, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FVEC(M), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. External Subroutines ..
      EXTERNAL          LSFUN1
C     .. Executable Statements ..
      CALL LSFUN1(M,N,X,FVEC)
      RETURN
C
C     END OF E04FDY   (GTVC)
C
      END
