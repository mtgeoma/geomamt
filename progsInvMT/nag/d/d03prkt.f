      SUBROUTINE D03PRK(XOP,NIP,XNP,NOP,ID,DXMESH)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C  ---------------------------------------------------------------------
C     METEST routine from SPRINT
C   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DXMESH
      INTEGER           ID, NIP, NOP
C     .. Array Arguments ..
      DOUBLE PRECISION  XNP(*), XOP(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, D2, D3, D4
      INTEGER           I, JJ, JK
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
      ID = 0
      DO 20 I = 1, NIP
         JJ = MIN(NIP,I+1)
         JK = MAX(1,I-1)
         D1 = XOP(JJ) - XOP(I)
         D2 = XOP(I) - XOP(JK)
         D3 = MAX(D1,D2)*DXMESH
         D4 = ABS(XOP(I)-XNP(I))
         IF (D4.GT.D3) GO TO 40
   20 CONTINUE
      RETURN
   40 ID = 1
      RETURN
      END
