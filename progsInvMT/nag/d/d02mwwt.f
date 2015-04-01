      SUBROUTINE D02MWW(N,H,HOLD,K,NMETH1,Y,DY,NY,W,THETA,IND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HOLD, THETA
      INTEGER           IND, K, N, NMETH1, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  DY(N), W(NY,4), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIF, FAC
      INTEGER           J
C     .. Executable Statements ..
C***********************************************************************
C  ROUTINE TO PREDICT THE NEW SOLUTION AND ITS TIME DERIVATIVE.
C***********************************************************************
      FAC = H/HOLD
      IF (K.EQ.1) THEN
C         INITIAL STEP OR RESTART FROM ORDER 1
         DO 20 J = 1, N
            DY(J) = W(J,2)/HOLD
            Y(J) = W(J,1) + DY(J)*H
   20    CONTINUE
      ELSE IF (NMETH1.EQ.30) THEN
C          NON-INITIAL STEP, PREDICT THE SOLUTION TO 2ND ORDER ACCURACY
C          USING INFORMATION FROM NEWTON ITERATION.
         DIF = 1.0D0 + THETA*(FAC-1.0D0)
         DO 40 J = 1, N
            Y(J) = W(J,1) + (DIF*W(J,3)+W(J,4))*FAC
            DY(J) = (W(J,2)+(DIF*W(J,3)+W(J,4)-W(J,2))/THETA)/HOLD
   40    CONTINUE
      ELSE
C          NON INITIAL STEP USE INFERIOR PREDICTOR BASED ON TRAPEZOIDAL
C          RULE AS FUNCTIONAL ITERATION IS BEING USED.
         DO 60 J = 1, N
            Y(J) = W(J,1) + FAC*W(J,2) + FAC**2*(W(J,2)-W(J,4))
            DY(J) = W(J,2)/HOLD + FAC/(THETA*HOLD)*(W(J,2)-W(J,4))
   60    CONTINUE
      END IF
      IF (IND.EQ.1) RETURN
      DO 80 J = 1, N
         W(J,1) = Y(J)
         W(J,2) = DY(J)*H
   80 CONTINUE
      RETURN
      END
