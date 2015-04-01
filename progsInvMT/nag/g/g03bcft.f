      SUBROUTINE G03BCF(STAND,PSCALE,N,M,X,LDX,Y,LDY,YHAT,R,LDR,ALPHA,
     *                  RSS,RES,WK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1112 (JUL 1993).
C
C     Computes Procrustes Rotations ( returned in R) for target matrix Y
C     for data in X and fitted values in YHAT
C
C     STAND indicates initial scaling/translation
C           STAND = 'N' no translation or scaling
C           STAND = 'Z' translation to origin
C           STAND = 'C' translation to Y centre
C           STAND = 'U' unit scaling, no translation
C           STAND = 'M' matched scaling and translation to Y
C           STAND = 'S' unit scaling and translation to zero
C
C     PSCALE indicates if post-rotation scaling is required
C            (Unscaled/Scaled)
C
C     ALPHA is the least-squares scaling coefficient
C
C     RSS is the residual sum of squares
C
C     RES is the distance between observed and fitted points
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03BCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, RSS
      INTEGER           IFAIL, LDR, LDX, LDY, M, N
      CHARACTER         PSCALE, STAND
C     .. Array Arguments ..
      DOUBLE PRECISION  R(LDR,M), RES(N), WK(M*M+7*M), X(LDX,M),
     *                  Y(LDY,M), YHAT(LDY,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  RN, SCALEF, SUM, SXX, SXY
      INTEGER           I, IERROR, IFAULT, J, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP1(1,1), WKSP2(1,1)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  F06RAF
      INTEGER           P01ABF
      EXTERNAL          F06RAF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02WEF, F06QFF, DGEMM, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (N.LE.0) THEN
         WRITE (P01REC,FMT=99999) N
      ELSE IF (M.LE.0) THEN
         WRITE (P01REC,FMT=99998) M
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC,FMT=99997) LDX, N
      ELSE IF (LDY.LT.N) THEN
         WRITE (P01REC,FMT=99996) LDY, N
      ELSE IF (LDR.LT.M) THEN
         WRITE (P01REC,FMT=99995) LDR, M
      ELSE IF (STAND.NE.'N' .AND. STAND.NE.'Z' .AND. STAND.NE.'C' .AND.
     *         STAND.NE.'U' .AND. STAND.NE.'S' .AND. STAND.NE.'M' .AND.
     *         STAND.NE.'n' .AND. STAND.NE.'z' .AND. STAND.NE.'c' .AND.
     *         STAND.NE.'u' .AND. STAND.NE.'s' .AND. STAND.NE.'m') THEN
         WRITE (P01REC,FMT=99994) STAND
      ELSE IF (PSCALE.NE.'S' .AND. PSCALE.NE.'U' .AND. PSCALE.NE.
     *         's' .AND. PSCALE.NE.'u') THEN
         WRITE (P01REC,FMT=99993) PSCALE
      ELSE
         IERROR = 0
      END IF
C
C     Translate and standardise if requested
C
      IF (IERROR.EQ.0) THEN
         IF (STAND.NE.'N' .AND. STAND.NE.'n' .AND. STAND.NE.'U' .AND.
     *       STAND.NE.'u') THEN
            RN = DBLE(N)
            DO 100 J = 1, M
               SUM = 0.0D0
               DO 20 I = 1, N
                  SUM = SUM + X(I,J)
   20          CONTINUE
               SUM = SUM/RN
               DO 40 I = 1, N
                  X(I,J) = X(I,J) - SUM
   40          CONTINUE
               SUM = 0.0D0
               DO 60 I = 1, N
                  SUM = SUM + Y(I,J)
   60          CONTINUE
               SUM = SUM/RN
               WK(J) = SUM
               DO 80 I = 1, N
                  Y(I,J) = Y(I,J) - SUM
   80          CONTINUE
  100       CONTINUE
         END IF
         IF (STAND.EQ.'U' .OR. STAND.EQ.'u' .OR. STAND.EQ.'S' .OR.
     *       STAND.EQ.'s' .OR. STAND.EQ.'M' .OR. STAND.EQ.'m') THEN
            SUM = F06RAF('F',N,M,X,LDX,WKSP1)
            IF (SUM.LE.0.0D0) THEN
               IERROR = 2
               WRITE (P01REC,FMT=99991)
               GO TO 300
            END IF
            SCALEF = 1.0D0/SUM
            SUM = F06RAF('F',N,M,Y,LDY,WKSP1)
            IF (SUM.LE.0.0D0) THEN
               IERROR = 2
               WRITE (P01REC,FMT=99990)
               GO TO 300
            END IF
            IF (STAND.EQ.'M' .OR. STAND.EQ.'m') THEN
               SCALEF = SUM*SCALEF
            ELSE
               SUM = 1.0D0/SUM
               DO 120 J = 1, M
                  CALL DSCAL(N,SUM,Y(1,J),1)
  120          CONTINUE
            END IF
            DO 140 J = 1, M
               CALL DSCAL(N,SCALEF,X(1,J),1)
  140       CONTINUE
         END IF
C
C        Compute SVD of X'Y = UDV'
C
         CALL DGEMM('T','N',M,M,N,1.0D0,X,LDX,Y,LDY,0.0D0,YHAT,LDY)
         IFAULT = -1
         CALL F02WEF(M,M,YHAT,LDY,0,WKSP1,1,.TRUE.,WKSP2,1,WK(M+1),
     *               .TRUE.,R,LDR,WK(2*M+1),IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 4
            WRITE (P01REC,FMT=99992)
            GO TO 300
         END IF
C
C        Compute rotation matrix R = UV' and rotated points YHAT = QR
C
         CALL F06QFF('G',M,M,R,LDR,WK(2*M+1),M)
         CALL DGEMM('N','N',M,M,M,1.0D0,YHAT,LDY,WK(2*M+1),M,0.0D0,R,
     *              LDR)
         CALL DGEMM('N','N',N,M,M,1.0D0,X,LDX,R,LDR,0.0D0,YHAT,LDY)
C
C        Compute scaling factor alpha
C
         IF (PSCALE.EQ.'S' .OR. PSCALE.EQ.'s') THEN
            SXX = F06RAF('F',N,M,YHAT,LDY,WKSP1)
            SXY = 0.0D0
            IF (SXX.LE.0.0D0) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99989)
               GO TO 300
            END IF
            DO 160 J = 1, M
               SXY = SXY + WK(M+J)
  160       CONTINUE
            ALPHA = SXY/(SXX*SXX)
            DO 180 J = 1, M
               CALL DSCAL(N,ALPHA,YHAT(1,J),1)
  180       CONTINUE
         END IF
C
C        Translate to original position if requested
C
         IF (STAND.EQ.'C' .OR. STAND.EQ.'c' .OR. STAND.EQ.'M' .OR.
     *       STAND.EQ.'m') THEN
            DO 240 J = 1, M
               SUM = WK(J)
               DO 200 I = 1, N
                  Y(I,J) = Y(I,J) + SUM
  200          CONTINUE
               DO 220 I = 1, N
                  YHAT(I,J) = YHAT(I,J) + SUM
  220          CONTINUE
  240       CONTINUE
         END IF
C
C        Compute residuals and RSS
C
         RSS = 0.0D0
         DO 280 I = 1, N
            SUM = 0.0D0
            DO 260 J = 1, M
               SUM = SUM + (Y(I,J)-YHAT(I,J))*(Y(I,J)-YHAT(I,J))
  260       CONTINUE
            RES(I) = SQRT(SUM)
            RSS = RSS + SUM
  280    CONTINUE
      END IF
  300 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.le.0: N = ',I16)
99998 FORMAT (' ** On entry, M.le.0: M = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N: LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, LDY.lt.N: LDY = ',I16,' N = ',I16)
99995 FORMAT (' ** On entry, LDR.lt.M: LDR = ',I16,' M = ',I16)
99994 FORMAT (' ** On entry, STAND not valid: STAND = ',A1)
99993 FORMAT (' ** On entry, PSCALE not valid: PSCALE = ',A1)
99992 FORMAT (' ** SVD failed to converge')
99991 FORMAT (' ** Only one distinct point (centred at zero) in X array'
     *       )
99990 FORMAT (' ** Only one distinct point (centred at zero) in Y array'
     *       )
99989 FORMAT (' ** Only zero points in YHAT array')
      END
