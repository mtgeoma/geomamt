      SUBROUTINE G05FSF(VK,N,T,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Subroutine generates random variates in the range
C     [-Pi,Pi] from a Von Mises distribution with
C     density proportional to exp(VK*COS(VMISES))
C     using best and fisher's method
C
C     VK=Parameter of distribution (0<REAL)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05FSF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  VK
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  T(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, C, F, PI, R, SS, VMISES, W
      INTEGER           I, IERROR, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF, X01AAF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, X01AAF, X02ALF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ACOS, COS, LOG, SQRT
C     .. Save statement ..
      SAVE              R, SS
C     .. Data statements ..
      DATA              SS/-1.0D0/
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
      IF (VK.LE.0.0D0) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) VK
      ELSE IF (N.LE.0) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) N
      END IF
      IF (IERROR.EQ.0) THEN
         IF (VK.NE.SS) THEN
            R = 1.0D0 + SQRT(1.0D0+4.0D0*VK*VK)
            R = (R-SQRT(R+R))/(VK+VK)
            IF (R.EQ.0.0D0) THEN
               R = 1.0D0/VK
            ELSE
               R = (1.0D0+R*R)/(R+R)
            END IF
            SS = VK
         END IF
         PI = X01AAF(0.0D0)
         BIG = X02ALF()
         DO 40 I = 1, N
   20       F = COS(PI*G05CAF(0.0D0))
            F = (1.0D0+R*F)/(R+F)
            IF (F.GT.1.0D0) THEN
               F = 1.0D0
            ELSE IF (F.LT.-1.0D0) THEN
               F = -1.0D0
            END IF
            C = VK*(R-F)
            IF (C.LE.0.0D0) GO TO 20
            W = G05CAF(0.0D0)
            IF (C*(2.0D0-C).LE.W) THEN
               IF (C.LE.BIG*W) THEN
                  IF (LOG(C/W)+1.0D0.LT.C) GO TO 20
               ELSE IF (W.GT.0.0D0) THEN
                  IF (LOG(C)-LOG(W)+1.0D0.LT.C) GO TO 20
               END IF
            END IF
            VMISES = ACOS(F)
            IF (G05CAF(0.0D0).LT.0.5D0) VMISES = -VMISES
            T(I) = VMISES
   40    CONTINUE
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, VK.le.0.0 : VK = ',D13.5)
99998 FORMAT (' ** On entry, N.le.0 : N = ',I16)
      END
