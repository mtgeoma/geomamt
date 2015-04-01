      SUBROUTINE G04CAT(ITAB,IDIM)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     AS 172.2
C
C     .. Scalar Arguments ..
      INTEGER           IDIM
C     .. Array Arguments ..
      INTEGER           ITAB(IDIM)
C     .. Local Scalars ..
      INTEGER           I, IK, ITEMP, ITER, K
C     .. Executable Statements ..
      ITER = IDIM/2
      K = IDIM + 1
      DO 20 I = 1, ITER
         ITEMP = ITAB(I)
         IK = K - I
         ITAB(I) = ITAB(IK)
         ITAB(IK) = ITEMP
   20 CONTINUE
      RETURN
      END
