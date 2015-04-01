      SUBROUTINE E04FDZ(M,N,X,FVEC,FJAC,LJ,S,IGRADE,NITER,NFTOTL,IW,LIW,
     *                  W,LW)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     DUMMY PRINT ROUTINE FOR THE EASY TO USE NONLINEAR LEAST-SQUARE
C     SUBROUTINES.
C
C     PHILIP E. GILL, SUSAN M. PICKEN AND WALTER MURRAY
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IGRADE, LIW, LJ, LW, M, N, NFTOTL, NITER
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), S(N), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Executable Statements ..
      RETURN
C
C     END OF E04FDZ   (LSMONA)
C
      END
