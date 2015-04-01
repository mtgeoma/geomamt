      SUBROUTINE G07DAF(N,X,Y,XME,XMD,XSD,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED. IER-692 (DEC 1989).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G07DAF')
      DOUBLE PRECISION  ZERO, PHI
      PARAMETER         (ZERO=0.0D0,PHI=0.6744897501962755D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMD, XME, XSD
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X1, X2
      INTEGER           I, IERROR, IFAIL2, K, K1, K2, KM, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      NREC = 1
C
C     PARAMETER CHECK
C
      IERROR = 0
      IF (N.EQ.2) THEN
         IF (X(2).GE.X(1)) THEN
            Y(1) = X(1)
            Y(2) = X(2)
         ELSE
            Y(1) = X(2)
            Y(2) = X(1)
         END IF
         XME = (Y(1)+Y(2))/2.D0
         XMD = Y(2) - XME
         XSD = XMD/PHI
      ELSE IF (N.GT.2) THEN
         DO 20 I = 1, N
            Y(I) = X(I)
   20    CONTINUE
C
C        Sort data.
C
         IFAIL2 = 0
         CALL M01CAF(Y,1,N,'A',IFAIL2)
         KM = (N+1)/2
         XME = Y(KM)
         IF (KM*2.EQ.N) XME = (XME+Y(KM+1))/2.D0
         K = 0
         K1 = KM
         K2 = KM
         X1 = ZERO
         X2 = ZERO
   40    CONTINUE
         IF (K.LT.KM) THEN
            K = K + 1
            IF (X1.GT.X2) THEN
               K2 = K2 + 1
               IF (K2.LE.N) THEN
                  X2 = Y(K2) - XME
                  GO TO 40
               END IF
            ELSE
               K1 = K1 - 1
               IF (K1.NE.0) THEN
                  X1 = XME - Y(K1)
                  GO TO 40
               END IF
            END IF
         END IF
         XMD = MIN(X1,X2)
         XSD = XMD/PHI
      ELSE
         IERROR = 1
         WRITE (REC,FMT=99999) N
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.le.1: N = ',I16)
      END
