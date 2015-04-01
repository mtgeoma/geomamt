      SUBROUTINE D02TKR(KD,MSTAR,N,V,Z,DMZ)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C**********************************************************************
C
C   purpose
C          compute dmz in a blockwise manner
C          dmz(i) = dmz(i)  +  v(i) * z(i), i = 1,...,n
C
C**********************************************************************
C
C
C     .. Scalar Arguments ..
      INTEGER           KD, MSTAR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DMZ(KD,N), V(KD,MSTAR,N), Z(MSTAR,N)
C     .. Local Scalars ..
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          DGEMV
C     .. Executable Statements ..
C
      DO 20 I = 1, N
         CALL DGEMV('n',KD,MSTAR,1.0D0,V(1,1,I),KD,Z(1,I),1,1.0D0,
     *              DMZ(1,I),1)
   20 CONTINUE
      RETURN
      END
