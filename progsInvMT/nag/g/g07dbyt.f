      DOUBLE PRECISION FUNCTION G07DBY(T,IPSI,DCHI)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Purpose
C     -------
C     Gives the value of the function CHI(T)=T*T/2 if PSI is 0 or
C     greater than 4, and CHI(T)=S*S/2 where S=MAX(-DCHI,MIN(DCHI,X))
C     otherwise.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DCHI, T
      INTEGER                          IPSI
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABST, PS
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MIN
C     .. Executable Statements ..
      IF (IPSI.GT.0 .AND. IPSI.LE.4) THEN
         ABST = ABS(T)
         PS = MIN(DCHI,ABST)
         G07DBY = PS*PS/2
      ELSE
         G07DBY = T*T/2
      END IF
      RETURN
      END
