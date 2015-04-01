      SUBROUTINE E04CGY(IFLAG,N,X,F,G,IW,LIW,W,LW)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     S-ROUTINE TO COMPUTE THE FUNCTION VALUE FOR E04JAF (I-ROUTINE
C     FOR E04JBL) AND ALSO FOR E04CGF.
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
C     .. Local Scalars ..
      DOUBLE PRECISION  H, XI
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          FUNCT1
C     .. Executable Statements ..
      IF (IFLAG.EQ.3) GO TO 20
      CALL FUNCT1(N,X,F)
      RETURN
   20 DO 40 I = 1, N
         H = G(I)
         XI = X(I)
         X(I) = XI + H
         CALL FUNCT1(N,X,G(I))
         X(I) = XI
   40 CONTINUE
      RETURN
C
C     END OF E04CGY (SBFUN1)
C
      END
