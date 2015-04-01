      SUBROUTINE G03EAF(UPDATE,DIST,SCALE,N,M,X,LDX,ISX,S,D,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes a distance matrix, D.
C
C     Three distances can be computed
C     DIST = A - absolute distances
C     DIST = E - Euclidean distances
C     DIST = S - Euclidean squared distance
C
C     SCALE indicates the scaling of the variables to be used
C     SCALE = S - standard deviation
C     SCALE = R - range
C     SCALE = G - scaling Given in S
C     SCALE = U - unscaled
C
C     The UPDATE option allows an existing matrix to be updated
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDX, M, N
      CHARACTER         DIST, SCALE, UPDATE
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N*(N-1)/2), S(M), X(LDX,M)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIJ, SUM, SW, XMAX, XMIN
      INTEGER           I, IERROR, IFAULT, IJ, J, K
      LOGICAL           EUCLID, ISXPOS, SRT
C     .. Local Arrays ..
      DOUBLE PRECISION  WMEAN(1), WT(1)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02BUF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
      IERROR = 0
      IF (N.LT.2) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) N
      ELSE IF (LDX.LT.N) THEN
         IERROR = 1
         WRITE (REC,FMT=99998) LDX, N
      ELSE IF (M.LE.0) THEN
         IERROR = 1
         WRITE (REC,FMT=99997) M
      ELSE IF (DIST.NE.'A' .AND. DIST.NE.'a' .AND. DIST.NE.'E' .AND.
     *         DIST.NE.'e' .AND. DIST.NE.'S' .AND. DIST.NE.'s') THEN
         IERROR = 1
         WRITE (REC,FMT=99996) DIST
      ELSE IF (UPDATE.NE.'I' .AND. UPDATE.NE.'i' .AND. UPDATE.NE.
     *         'U' .AND. UPDATE.NE.'u') THEN
         IERROR = 1
         WRITE (REC,FMT=99995) UPDATE
      ELSE IF (SCALE.NE.'S' .AND. SCALE.NE.'s' .AND. SCALE.NE.'R' .AND.
     *         SCALE.NE.'r' .AND. SCALE.NE.'G' .AND. SCALE.NE.'g' .AND.
     *         SCALE.NE.'U' .AND. SCALE.NE.'u') THEN
         IERROR = 1
         WRITE (REC,FMT=99994) SCALE
      END IF
      IF (IERROR.EQ.0) THEN
         ISXPOS = .FALSE.
         DO 20 I = 1, M
            IF (ISX(I).GT.0) ISXPOS = .TRUE.
   20    CONTINUE
         IF ( .NOT. ISXPOS) THEN
            IERROR = 2
            WRITE (REC,FMT=99993)
         ELSE
            IF (UPDATE.EQ.'I' .OR. UPDATE.EQ.'i') THEN
               DO 40 I = 1, N*(N-1)/2
                  D(I) = 0.0D0
   40          CONTINUE
            ELSE
               DO 60 I = 1, N*(N-1)/2
                  IF (D(I).LT.0.0D0) THEN
                     IERROR = 2
                     WRITE (REC,FMT=99992)
                     GO TO 240
                  END IF
   60          CONTINUE
            END IF
            IF (DIST.EQ.'A' .OR. DIST.EQ.'a') THEN
               EUCLID = .FALSE.
               SRT = .FALSE.
            ELSE IF (DIST.EQ.'E' .OR. DIST.EQ.'e') THEN
               EUCLID = .TRUE.
               SRT = .TRUE.
            ELSE
               EUCLID = .TRUE.
               SRT = .FALSE.
            END IF
            IF (SCALE.EQ.'S' .OR. SCALE.EQ.'s') THEN
               DO 80 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     IFAULT = -1
                     CALL G02BUF('MEAN','UNWEIGHTED',N,1,X(1,J),LDX,WT,
     *                           SW,WMEAN,S(J),IFAULT)
                     IF (S(J).LE.0.0D0) THEN
                        IERROR = 2
                        WRITE (REC,FMT=99991) J
                        GO TO 240
                     END IF
                     S(J) = SQRT(S(J)/(SW-1.0D0))
                  END IF
   80          CONTINUE
            ELSE IF (SCALE.EQ.'R' .OR. SCALE.EQ.'r') THEN
               DO 120 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     XMAX = X(1,J)
                     XMIN = X(1,J)
                     DO 100 I = 2, N
                        IF (X(I,J).GT.XMAX) XMAX = X(I,J)
                        IF (X(I,J).LT.XMIN) XMIN = X(I,J)
  100                CONTINUE
                     IF (XMIN.GE.XMAX) THEN
                        IERROR = 2
                        WRITE (REC,FMT=99991) J
                        GO TO 240
                     END IF
                     S(J) = XMAX - XMIN
                  END IF
  120          CONTINUE
            ELSE IF (SCALE.EQ.'G' .OR. SCALE.EQ.'g') THEN
               DO 140 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     IF (S(J).LE.0.0D0) THEN
                        IERROR = 2
                        WRITE (REC,FMT=99990)
                        GO TO 240
                     END IF
                  END IF
  140          CONTINUE
            ELSE
               DO 160 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     S(J) = 1.0D0
                  END IF
  160          CONTINUE
            END IF
            IJ = 0
            DO 220 I = 1, N
               DO 200 J = 1, I - 1
                  IJ = IJ + 1
                  SUM = 0.0D0
                  DO 180 K = 1, M
                     IF (ISX(K).GT.0) THEN
                        DIJ = (X(I,K)-X(J,K))/S(K)
                        IF (EUCLID) THEN
                           SUM = SUM + DIJ*DIJ
                        ELSE
                           SUM = SUM + ABS(DIJ)
                        END IF
                     END IF
  180             CONTINUE
                  IF (SRT) SUM = SQRT(SUM)
                  D(IJ) = D(IJ) + SUM
  200          CONTINUE
  220       CONTINUE
         END IF
      END IF
  240 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.2: N = ',I16)
99998 FORMAT (1X,'** On entry, LDX.lt.N: LDX = ',I16,' N = ',I16)
99997 FORMAT (1X,'** On entry, M.le.0: M = ',I16)
99996 FORMAT (1X,'** On entry, DIST is not valid: DIST = ',A1)
99995 FORMAT (1X,'** On entry, UPDATE is not valid: UPDATE = ',A1)
99994 FORMAT (1X,'** On entry, SCALE is not valid: SCALE = ',A1)
99993 FORMAT (1X,'** On entry, ISX does not contain a positive element')
99992 FORMAT (1X,'** On entry, at least one element of D.lt.0.0')
99991 FORMAT (1X,'** Variable ',I16,' is constant.')
99990 FORMAT (1X,'** On entry, at least one element of S.le.0.0')
      END
