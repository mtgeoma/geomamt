      INTEGER FUNCTION G01ARS(EPSI,MAXINT,X,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-917 (APR 1991).
C
C     Find the integer equal to or next closer to zero than X.
C
C     Tests whether  X  is too large to fit in an integer variable.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        EPSI, MAXINT, X
      INTEGER                 IFAIL
C     .. External Functions ..
      INTEGER                 X02BBF
      EXTERNAL                X02BBF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, INT
C     .. Executable Statements ..
C
      IF (ABS(X).LT.MAXINT) THEN
C
         G01ARS = INT((1.0D0+10.0D0*EPSI)*X)
      ELSE
C
C        X  is too large in magnitude to fit in an integer.
C        Return the largest legal integer and set the error flag.
C
         IFAIL = 4
         G01ARS = X02BBF(0.0D0)
      END IF
      RETURN
      END
