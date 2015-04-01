      SUBROUTINE E04EBZ(IFLAG,N,X,HESL,LH,HESD,IW,LIW,W,LW)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     S-ROUTINE TO COMPUTE HESSIAN MATRIX FOR E04LAF (I-ROUTINE FOR
C     E04LBR) AND ALSO FOR E04EBF.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IFLAG, LH, LIW, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  HESD(N), HESL(LH), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. External Subroutines ..
      EXTERNAL          HESS2
C     .. Executable Statements ..
      CALL HESS2(N,X,HESL,LH,HESD)
      RETURN
C
C     END OF E04EBZ (SBHESS)
C
      END
