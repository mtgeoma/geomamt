      SUBROUTINE D02NMR(ABD,LDA,N,ML,MU,RDAE,H,EL0)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     OLD NAME BANDSC
C
C***********************************************************************
C     ROUTINE TO SCALE THE ROWS OF THE MATRIX IN ABD BY 1/(H*EL0)
C     WHERE H*EL0 IS A FACTOR
C***********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EL0, H
      INTEGER           LDA, ML, MU, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ABD(LDA,N), RDAE(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EL0H, SCALE
      INTEGER           I, J, M
C     .. Executable Statements ..
      EL0H = 1.0D0/(EL0*H)
      M = ML + MU + 1
      DO 40 I = 1, N
         SCALE = RDAE(I) + (1.0D0-RDAE(I))*EL0H
         DO 20 J = 1, LDA
            ABD(J,I) = SCALE*ABD(J,I)
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
