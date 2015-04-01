      SUBROUTINE D02QDW(X,XEND,Y,N,IW,CIN,COMM,COUT,ISAVE,CHK1,CHK2)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, XEND
      INTEGER           ISAVE, IW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CHK1(N), CHK2(N), CIN(7), COMM(5), COUT(16),
     *                  Y(N)
C     .. Scalars in Common ..
      INTEGER           IWCOPY, NCOPY
C     .. Local Scalars ..
      DOUBLE PRECISION  DIR, DIRCM5
      INTEGER           I
      LOGICAL           ORDROK
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN
C     .. Common blocks ..
      COMMON            /AD02QQ/NCOPY, IWCOPY
C     .. Save statement ..
      SAVE              /AD02QQ/
C     .. Executable Statements ..
      IF (COMM(1).LT.0.0D0) THEN
         CIN(1) = -9.0D0
         GO TO 60
      END IF
C
      IF (COMM(2).LT.0.0D0) THEN
         CIN(1) = -10.0D0
         GO TO 60
      ELSE IF (COMM(2).GT.0.0D0) THEN
         DO 20 I = 1, N
            IF (ABS(Y(I)).LT.COMM(2)) GO TO 20
            CIN(1) = -11.0D0
            COMM(2) = 0.0D0
            GO TO 60
   20    CONTINUE
      END IF
C
      IF (COMM(3).NE.0.0D0) THEN
         IF (N.NE.NCOPY .OR. IW.NE.IWCOPY) THEN
            CIN(1) = -19.0D0
            GO TO 60
         END IF
         DO 40 I = 1, N
            CHK2(I) = Y(I)
            IF (Y(I).NE.CHK1(I)) GO TO 40
            CIN(1) = -12.0D0
            COMM(3) = 0.0D0
            GO TO 60
   40    CONTINUE
      END IF
C
      IF (COMM(4).LT.0.0D0) THEN
         DIR = SIGN(1.D0,XEND-X)
         DIRCM5 = DIR*COMM(5)
         IF (CIN(1).EQ.1.0D0) THEN
            ORDROK = DIR*X .LE. DIRCM5 .AND. DIRCM5 .LE. XEND
            IF ( .NOT. ORDROK) THEN
               CIN(1) = -13.0D0
               COMM(4) = 0.0D0
            END IF
         ELSE IF (CIN(1).GT.1.0D0) THEN
            ORDROK = DIR*COUT(4) .LT. DIR*X .AND. DIR*X .LT. DIR*XEND
            IF (COUT(8).GT.1.0D0) ORDROK = DIR*COUT(5) .LT. DIR*COUT(4)
     *           .AND. ORDROK
            IF ( .NOT. ORDROK) THEN
               CIN(1) = -14.0D0
               COMM(4) = 0.0D0
               GO TO 60
            END IF
            ORDROK = DIR*COUT(4) .LE. DIRCM5 .AND. DIRCM5 .LE. DIR*XEND
            IF ( .NOT. ORDROK) THEN
               CIN(1) = -15.0D0
               COMM(4) = 0.0D0
            ELSE IF (DIRCM5.LE.DIR*X) THEN
               ISAVE = 0
               CIN(1) = 6.0D0
               COMM(4) = 0.0D0
            END IF
         END IF
      END IF
   60 CONTINUE
      RETURN
      END
