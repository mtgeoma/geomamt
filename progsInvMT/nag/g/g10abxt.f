      SUBROUTINE G10ABX(X,AVH,Y,DY,N,Q,YHAT,C,LDC,SU,SV)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Calculates coefficients of a cubic smoothing spline from
C     parameters calculated by subroutine G10ABY.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AVH, Q
      INTEGER           LDC, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,3), DY(N), SU(N+2), SV(N+2), X(N), Y(N),
     *                  YHAT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, QH
      INTEGER           I
C     .. Executable Statements ..
C
C     Calculate A
C
      QH = Q/(AVH*AVH)
      DO 20 I = 1, N
         SU(I+1) = QH*SU(I+1)
   20 CONTINUE
C
C     Calculate C
C
      DO 40 I = 1, N - 1
         H = X(I+1) - X(I)
         C(I,3) = (SU(I+2)-SU(I+1))/(3.0D0*H)
         C(I,1) = (YHAT(I+1)-YHAT(I))/H - (H*C(I,3)+SU(I+1))*H
         C(I,2) = SU(I+1)
   40 CONTINUE
      RETURN
      END
