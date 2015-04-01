      SUBROUTINE E04HCZ(N,Y,Z)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     FORMS TWO ORTHOGONAL VECTORS Y AND Z. THE ELEMENTS OF BOTH
C     VECTORS ARE LESS THAN 1.0 IN ABSOLUTE VALUE. THE ELEMENTS OF Z
C     ARE POSITIVE AND Z(1) = 0.5726. IF N .EQ. 1 Y(1) = 0.0.
C
C     PHILIP E. GILL, ENID M. R. LONG, WALTER MURRAY AND SUSAN M.
C     PICKEN, D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N), Z(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, REVALP, WI, ZI, ZSUM, ZWSUM
      INTEGER           I
C     .. Executable Statements ..
      ZI = 0.5726D+0
      Z(1) = ZI
      Y(1) = (1.0D+0-ZI)**2
      ZSUM = ZI*ZI
      ZWSUM = ZI*Y(1)
      IF (N.GE.2) GO TO 20
      Y(1) = 0.0D+0
      RETURN
C
C     FORM Z(I) = Z(I - 1)**2 AND ADD 0.88 IF LESS THAN 0.1. FORM
C     W(I) = (1.0 - Z(I))**2 AND STORE IN Y. FORM ALPHA = ZT*Z/ZT*W.
C
   20 DO 40 I = 2, N
         ZI = ZI*ZI
         IF (ZI.LT.0.1D+0) ZI = ZI + 0.88D+0
         WI = (1.0D+0-ZI)**2
         Y(I) = WI
         Z(I) = ZI
         ZSUM = ZSUM + ZI*ZI
         ZWSUM = ZWSUM + ZI*WI
   40 CONTINUE
      ALPHA = ZSUM/ZWSUM
      REVALP = 1.0D+0/(1.0D+0+ALPHA)
C
C     FORM Y(I) = (Z(I) - ALPHA*W(I))/(1.0 + ALPHA).
C
      DO 60 I = 1, N
         Y(I) = REVALP*(Z(I)-ALPHA*Y(I))
   60 CONTINUE
      RETURN
C
C     END OF E04HCZ (ORTHVC)
C
      END
