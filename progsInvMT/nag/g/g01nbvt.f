      DOUBLE PRECISION FUNCTION G01NBV(X,IWK,WK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C      Computes function for integration by G01NBT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Array Arguments ..
      DOUBLE PRECISION                 WK(*)
      INTEGER                          IWK(*)
C     .. Local Scalars ..
      INTEGER                          ICODE, IMU, IS, ISROW, LWKAA,
     *                                 LWKCC, LWKR, LWKSLA, LWKSMU, N,
     *                                 NN
C     .. External Functions ..
      DOUBLE PRECISION                 G01NBZ
      EXTERNAL                         G01NBZ
C     .. Executable Statements ..
      ICODE = IWK(1)
      IMU = IWK(2)
      N = IWK(3)
      NN = N*(N+1)/2
      IS = IWK(4)
      ISROW = IWK(5)
      LWKR = IWK(6)
      LWKSMU = IWK(7)
      LWKSLA = IWK(8)
      LWKAA = IWK(9)
      LWKCC = IWK(10)
      G01NBV = G01NBZ(X,ICODE,IMU,N,IS,ISROW,IWK(11),WK(LWKR),WK(LWKSMU)
     *         ,WK(LWKSLA),WK(LWKAA),WK(LWKCC),WK,WK(N+1),WK(2*N+1),
     *         WK(2*N+NN+1),WK(2*N+2*NN+1),WK(2*N+3*NN+1),WK(2*N+4*NN+1)
     *         )
      RETURN
      END
