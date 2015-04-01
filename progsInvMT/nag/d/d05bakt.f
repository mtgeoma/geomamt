      SUBROUTINE D05BAK(CK,H,THETA,VK,WVK,IQ,NSTART)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     -----------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++
C       this routine evaluates the kernel associated with
C       a Runge-Kutta scheme.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++
C     -----------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H
      INTEGER           IQ, NSTART
C     .. Array Arguments ..
      DOUBLE PRECISION  THETA(0:IQ-1), VK(0:NSTART),
     *                  WVK(0:NSTART-1,0:IQ-1,0:IQ-1)
C     .. Function Arguments ..
      DOUBLE PRECISION  CK
      EXTERNAL          CK
C     .. Local Scalars ..
      INTEGER           I, IQM1, IR, IS, J
C     .. Executable Statements ..
C
      IQM1 = IQ - 1
      DO 20 I = 0, NSTART - 1
         WVK(I,0,0) = VK(I)
   20 CONTINUE
C
      DO 40 J = 1, IQM1
         WVK(0,J,J) = VK(0)
   40 CONTINUE
      DO 80 J = 1, IQM1
         DO 60 I = 1, NSTART - 1
            WVK(I,J,J) = VK(I)
            WVK(I,J,0) = CK(I*H+THETA(J)*H)
            WVK(I,0,J) = CK(I*H-THETA(J)*H)
   60    CONTINUE
   80 CONTINUE
C
      DO 120 IR = 1, IQM1
         DO 100 IS = 0, IR - 1
            WVK(0,IR,IS) = CK((THETA(IR)-THETA(IS))*H)
  100    CONTINUE
  120 CONTINUE
C
      DO 180 IS = 1, IQM1
         DO 160 IR = 1, IQM1
            IF (IR-IS.NE.0) THEN
               DO 140 I = 1, NSTART - 1
                  WVK(I,IR,IS) = CK(I*H+(THETA(IR)-THETA(IS))*H)
  140          CONTINUE
            END IF
  160    CONTINUE
  180 CONTINUE
C
      RETURN
      END
