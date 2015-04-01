      SUBROUTINE D02NMT(A,LDA,N,RDAE,H,EL0)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     OLD NAME FULLSC
C
C***********************************************************************
C     ROUTINE TO SCALE THE FULL MATRIX A(LDA,N) BY THE FACTOR
C     (H*EL0)**(-1) WHERE POSSIBLE
C***********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EL0, H
      INTEGER           LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), RDAE(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EL0H
      INTEGER           I, J
C     .. Executable Statements ..
      EL0H = 1.0D0/(EL0*H)
      DO 40 J = 1, N
         DO 20 I = 1, N
            A(I,J) = A(I,J)*(RDAE(I)+(1.0D0-RDAE(I))*EL0H)
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
