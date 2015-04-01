      SUBROUTINE D05BAY(A,B,THETA,IQ)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      INTEGER           IQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(6,0:5), B(0:6), THETA(0:6)
C     .. Executable Statements ..
C
      IQ = 7
C
      B(0) = 35.D0/384.D0
      B(1) = 0.0D0
      B(2) = 500.D0/1113.D0
      B(3) = 125.D0/192.D0
      B(4) = -2187.D0/6784.D0
      B(5) = 11.D0/84.D0
      B(6) = 0.D0
C
      THETA(0) = 0.0D0
      THETA(1) = 0.2D0
      THETA(2) = 3.D0/10.D0
      THETA(3) = 0.8D0
      THETA(4) = 8.D0/9.D0
      THETA(5) = 1.D0
      THETA(6) = 1.D0
C
      A(1,0) = 0.2D0
      A(2,0) = 3.D0/40.D0
      A(2,1) = 9.D0/40.D0
      A(3,0) = 44.D0/45.D0
      A(3,1) = -56.D0/15.D0
      A(3,2) = 32.D0/9.D0
      A(4,0) = 19372.D0/6561.D0
      A(4,1) = -25360.D0/2187.D0
      A(4,2) = 64448.D0/6561.D0
      A(4,3) = -212.D0/729.D0
      A(5,0) = 9017.D0/3168.D0
      A(5,1) = -355.D0/33.D0
      A(5,2) = 46732.D0/5247.D0
      A(5,3) = 49.D0/176.D0
      A(5,4) = -5103.D0/18656.D0
      A(6,0) = B(0)
      A(6,1) = B(1)
      A(6,2) = B(2)
      A(6,3) = B(3)
      A(6,4) = B(4)
      A(6,5) = B(5)
C
      RETURN
      END
