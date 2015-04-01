      SUBROUTINE E04DEZ(IFLAG,N,X,F,G,IW,LIW,W,LW)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     S-ROUTINE TO COMPUTE FUNCTION VALUE AND GRADIENT VECTOR FOR
C     E04LAF, E04KCF AND E04KAF (I-ROUTINES FOR E04LBR AND E04KBN)
C     AND ALSO FOR E04EBF, E04DFF AND E04DEF.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  F
      INTEGER           IFLAG, LIW, LW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. External Subroutines ..
      EXTERNAL          FUNCT2
C     .. Executable Statements ..
      CALL FUNCT2(N,X,F,G)
      RETURN
C
C     END OF E04DEZ (SBFUN2)
C
      END
