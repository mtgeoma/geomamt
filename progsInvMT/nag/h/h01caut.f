      SUBROUTINE H01CAU(N,V)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     H01CAU NEGATIVES EVERY ELEMENT OF THE VECTOR V.
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY AND
C     SUSAN M. PICKEN
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND
C     JULY 1977
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  V(N)
C     .. Local Scalars ..
      INTEGER           I
C     .. Executable Statements ..
      DO 20 I = 1, N
         V(I) = -V(I)
   20 CONTINUE
      RETURN
C
C     END OF H01CAU (NEGVEC)
C
      END
