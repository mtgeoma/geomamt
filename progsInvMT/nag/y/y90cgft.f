      SUBROUTINE Y90CGF(DET,VTYPE,N,V,VBOUND,COND,SCALE,DETMAN,DETEXP,
     *                  DIST,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ========================================
C         *  Y90CGF :  Complex Vector Generator  *
C         ========================================
C
C
C              ARGUMENTS
C              =========
C
C     DET    :  'D'  ==> Calculate the product of all the
C                        elements of the vector V.
C               any  ==> Do not calculate.
C
C     VTYPE  :  Defines how the vector V is generated.
C               +/-1 ==> the elements of V have values with real parts
C                        between the real parts of VBOUND(1),VBOUND(2),
C                        likewise for the imaginary parts
C               +/-2 ==> use the parameter COND to generate real
C                        elements using  V(I) = COND**(-(I-1)/(N-1))
C               +/-3 ==> as above using (1-(I-1)/(N+1)+(I-1)/(COND*(N-1)
C               +/-4 ==> as for +/-1 but real value are multiplied by
C                        random complex numbers on the unit circle
C               +/-5 ==> as for +/-2 but real value are multiplied by
C                        random complex numbers on the unit circle
C               If VTYPE < 0 the order of the entries of V is reversed.
C
C     N      :  Number of elements of the vector V
C
C     V      :  Generated vector V
C
C     VBOUND :  Contains on entry limiting values for the elements of V.
C
C     SCALE  :  Scaling factor which multiplies all elements of V.
C
C     DETMAN :  To avoid overflow the product of the elements of V
C               is calculated as: DETMAN*2.0**DETEXP, such that
C               0.0625 <= |DETMAN| < 1.0.
C
C     DETEXP :  Exponent used to calculate the product of the
C               elements of V (see DETMAN above).
C
C     DIST   :  Type of distribution used to generate random numbers.
C               For real random numbers (|VTYPE| = 1):
C               1 ==>  Uniform (0,1)
C               2 ==>  Uniform (-1,1)
C               3 ==>  Normal (0,1)
C               For complex random numbers (|VTYPE| = 4 or 5):
C               1 ==>  real and imaginary parts uniform (0,1)
C               2 ==>      "        "          "        (-1,1)
C               3 ==>  random comple in the disc |z| < 1
C               4 ==>  real and imaginary parts normal (0,1)
C               5 ==>  random complex on the circle |z| = 1
C
C     SEED   :  Seeds for the random number generator.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF
      COMPLEX*16        CZERO, CONE
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,CZERO=(0.0D0,0.0D0),
     *                  CONE=(1.0D0,0.0D0))
C     .. Scalar Arguments ..
      COMPLEX*16        DETMAN, SCALE
      DOUBLE PRECISION  COND
      INTEGER           DETEXP, DIST, N, VTYPE
      CHARACTER*1       DET
C     .. Array Arguments ..
      COMPLEX*16        V(*), VBOUND(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      DOUBLE PRECISION  XA, XB
      INTEGER           I
C     .. External Functions ..
      COMPLEX*16        Y90EBF
      DOUBLE PRECISION  Y90TBF
      LOGICAL           Y90WAF
      EXTERNAL          Y90EBF, Y90TBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          Y90DMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, DBLE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Generate vector
C
C-----------------------------------------------------------------------
      IF (ABS(VTYPE).EQ.1) THEN
         DO 20 I = 1, N
            IF (DIST.LE.1) THEN
               XA = DBLE(VBOUND(1)) + DBLE(VBOUND(2)-VBOUND(1))
     *              *Y90TBF(DIST,SEED)
               XB = DIMAG(VBOUND(1)) + DIMAG(VBOUND(2)-VBOUND(1))
     *              *Y90TBF(DIST,SEED)
            ELSE IF (DIST.EQ.2) THEN
               XA = HALF*(DBLE(VBOUND(1)+VBOUND(2))+DBLE(VBOUND(2)
     *              -VBOUND(1))*Y90TBF(DIST,SEED))
               XB = HALF*(DIMAG(VBOUND(1)+VBOUND(2))+DIMAG(VBOUND(2)
     *              -VBOUND(1))*Y90TBF(DIST,SEED))
            ELSE
               XA = DBLE(VBOUND(2))*Y90TBF(DIST,SEED) + DBLE(VBOUND(1))
               XB = DIMAG(VBOUND(2))*Y90TBF(DIST,SEED) + DIMAG(VBOUND(1)
     *              )
            END IF
            V(I) = DCMPLX(XA,XB)
   20    CONTINUE
      ELSE IF (ABS(VTYPE).EQ.2 .OR. ABS(VTYPE).EQ.4) THEN
         V(1) = CONE
         DO 40 I = 2, N
            V(I) = DCMPLX(COND**(-DBLE(I-1)/DBLE(N-1)),ZERO)
   40    CONTINUE
      ELSE IF (ABS(VTYPE).EQ.3 .OR. ABS(VTYPE).EQ.5) THEN
         V(1) = CONE
         DO 60 I = 2, N
            V(I) = DCMPLX(DBLE(1-DBLE(I-1)/DBLE(N-1))+(DBLE(I-1)
     *             /DBLE(N-1))/COND,ZERO)
   60    CONTINUE
      END IF
C
      IF (ABS(VTYPE).EQ.4 .OR. ABS(VTYPE).EQ.5) THEN
         DO 80 I = 1, N
            V(I) = V(I)*Y90EBF(5,SEED)
   80    CONTINUE
      END IF
C
      IF (ABS(VTYPE).NE.1) THEN
         DO 100 I = 1, N
            V(I) = V(I)*SCALE
  100    CONTINUE
      END IF
C
      IF (VTYPE.LE.-2) THEN
         DO 120 I = 1, N/2
            TEMP = V(I)
            V(I) = V(N-I+1)
            V(N-I+1) = TEMP
  120    CONTINUE
      END IF
C
C     Calculate the determinant
C
      IF (Y90WAF(DET,'D')) THEN
         DETMAN = CONE
         DETEXP = 0
         DO 140 I = 1, N
            DETMAN = DETMAN*V(I)
            CALL Y90DMF(DETMAN,DETEXP,4)
  140    CONTINUE
      ELSE
         DETMAN = CZERO
         DETEXP = 0
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90CGF
C
C-----------------------------------------------------------------------
      RETURN
      END
