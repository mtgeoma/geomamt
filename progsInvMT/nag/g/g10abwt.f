      SUBROUTINE G10ABW(X,AVH,N,R,P,H,W)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Calculates the diagonal elements of the influence matrix.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AVH, P
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  H(N), R(N+2,3), W(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, F1, G, G1, H1, ONE, S, ZERO
      INTEGER           I
C     .. Data statements ..
      DATA              ZERO, ONE/0.0D0, 1.0D0/
C     .. Executable Statements ..
C
C     Initialize
C
      S = AVH/(X(2)-X(1))
      H(1) = ONE - P*W(1)*W(1)*S*S*R(3,1)
      R(2,1) = ZERO
      R(2,2) = ZERO
      R(2,3) = ZERO
C
C     Calculate diagonal elements
C
      DO 20 I = 2, N - 1
         F = S
         S = AVH/(X(I+1)-X(I))
         G = -F - S
         F1 = F*R(I,1) + G*R(I,2) + S*R(I,3)
         G1 = F*R(I,2) + G*R(I+1,1) + S*R(I+1,2)
         H1 = F*R(I,3) + G*R(I+1,2) + S*R(I+2,1)
         H(I) = ONE - P*W(I)*W(I)*(F*F1+G*G1+S*H1)
   20 CONTINUE
      H(N) = ONE - P*W(N)*W(N)*S*S*R(N,1)
      RETURN
      END
