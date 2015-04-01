      DOUBLE PRECISION FUNCTION E02GBJ(NUM,VX,IX,VY,IY,KX,KY)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     *********************
C     INNER PRODUCT OF THE VECTORS VX AND VY.
C     THE COMPONENTS OF THE VECTOR VX ARE ASSUMED
C     TO BE STORED  AT EACH IX-TH POSITION ALONG THE ARRAY.
C     THE COMPONENTS OF VY ARE AT EACH IY-TH POSITION.
C     *********************
C
C     .. Scalar Arguments ..
      INTEGER                          IX, IY, KX, KY, NUM
C     .. Array Arguments ..
      DOUBLE PRECISION                 VX(KX), VY(KY)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SS, ZERO
      INTEGER                          JX, JXS, JY
C     .. Data statements ..
      DATA                             ZERO/0.0D+00/
C     .. Executable Statements ..
      JXS = ((NUM-1)*IX) + 1
      SS = ZERO
      IF (1.GT.NUM) GO TO 40
      JY = 1
      DO 20 JX = 1, JXS, IX
         SS = SS + (VX(JX))*(VY(JY))
         JY = JY + IY
   20 CONTINUE
   40 CONTINUE
      E02GBJ = SS
      RETURN
      END
