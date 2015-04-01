      SUBROUTINE D02MWV(N,H,HOLD,K,NMETH1,Y,DY,NY,W,THETA,IND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1541 (JUN 1995).
C     VP changed
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HOLD, THETA
      INTEGER           IND, K, N, NMETH1, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  DY(N), W(NY,4), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIF, FAC, FAC2
      INTEGER           J
C     .. Executable Statements ..
C   -------------------------------------------------------------------
C   | ROUTINE TO RETRACT THE SOLUTION AND ITS TIME DERIV TO THEIR     |
C   | VALUES  AT THE PREVIOUS TIME STEP.                              |
C   -------------------------------------------------------------------
      FAC = H/HOLD
      IF (K.EQ.1) THEN
         DO 20 J = 1, N
            W(J,1) = W(J,1) - W(J,2)
            W(J,2) = W(J,2)/FAC
   20    CONTINUE
      ELSE IF (NMETH1.EQ.30) THEN
         DIF = 1.0D0 + THETA*(FAC-1.0D0)
         FAC2 = -(1.0D0-THETA)*FAC
         DO 40 J = 1, N
C  VP CHANGED NEXT LINE SEPT 92
            W(J,1) = W(J,1) - (DIF*W(J,3)+W(J,4))*FAC
C           W(J,1) = W(J,1) - DIF*W(J,3) - W(J,4)*FAC
            W(J,2) = (W(J,2)*THETA-(DIF*W(J,3)+W(J,4))*FAC)/FAC2
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            W(J,2) = (W(J,2)+FAC**2/THETA*W(J,4))/(FAC+FAC**2/THETA)
            W(J,1) = W(J,1) - FAC*W(J,2) - FAC**2*(W(J,2)-W(J,4))
   60    CONTINUE
      END IF
      IF (IND.EQ.1) RETURN
      DO 80 J = 1, N
         Y(J) = W(J,1)
         DY(J) = W(J,2)/HOLD
   80 CONTINUE
      RETURN
      END
