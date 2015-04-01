      SUBROUTINE G02BWF(M,R,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES CORRELATION MATRIX FROM SUM OF
C     SQUARES AND CROSS-PRODUCTS OF DEVIATIONS
C     STORED IN PACKED FORM UPPER-TRIANGULAR BY COLUMNS
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BWF')
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  R((M*M+M)/2)
C     .. Local Scalars ..
      DOUBLE PRECISION  RS
      INTEGER           I, IC, IERROR, IJ, J, JC
      LOGICAL           NOWT
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IERROR = 1
      IF (M.LT.1) THEN
         WRITE (REC,FMT=99999) M
      ELSE
         NOWT = .FALSE.
         IERROR = 0
         IC = 0
         DO 20 I = 1, M
            IC = IC + I
            IF (R(IC).GT.ZERO) THEN
               R(IC) = ONE/SQRT(R(IC))
            ELSE
               NOWT = .TRUE.
            END IF
   20    CONTINUE
         IF (NOWT) THEN
            IERROR = 2
            WRITE (REC,FMT=99998)
            IJ = 0
            IC = 0
            DO 80 I = 1, M
               IC = IC + I
               RS = R(IC)
               IF (RS.GT.ZERO) THEN
                  JC = 0
                  DO 40 J = 1, I - 1
                     JC = JC + J
                     IJ = IJ + 1
                     R(IJ) = R(IJ)*RS*R(JC)
   40             CONTINUE
               ELSE
                  DO 60 J = 1, I - 1
                     IJ = IJ + 1
                     R(IJ) = ZERO
   60             CONTINUE
               END IF
               IJ = IJ + 1
   80       CONTINUE
            JC = 0
            DO 100 J = 1, M
               JC = JC + J
               IF (R(JC).GT.ZERO) R(JC) = ONE
  100       CONTINUE
         ELSE
            IJ = 0
            IC = 0
            DO 140 I = 1, M
               IC = IC + I
               RS = R(IC)
               JC = 0
               DO 120 J = 1, I - 1
                  IJ = IJ + 1
                  JC = JC + J
                  R(IJ) = R(IJ)*RS*R(JC)
  120          CONTINUE
               IJ = IJ + 1
  140       CONTINUE
            JC = 0
            DO 160 J = 1, M
               JC = JC + J
               R(JC) = ONE
  160       CONTINUE
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
C
99999 FORMAT (' ** On entry, M.lt.1 : M = ',I16)
99998 FORMAT (' ** On entry, a variable has zero ''variance''.')
      END
