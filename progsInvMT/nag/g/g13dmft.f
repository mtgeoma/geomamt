      SUBROUTINE G13DMF(MATRIX,K,N,M,W,IK,WMEAN,R0,R,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16A REVISED. IER-1044 (JUN 1993).
C     MARK 17 REVISED. IER-1691 (JUN 1995).
C
C     This subroutine calculates the sample cross-covariance
C     (MATRIX = 'V') or cross-correlation (MATRIX = 'R')
C     matrices of a multivariate time series.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      CHARACTER*6       SRNAME
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0,SRNAME='G13DMF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, K, M, N
      CHARACTER         MATRIX
C     .. Array Arguments ..
      DOUBLE PRECISION  R(IK,IK,M), R0(IK,K), W(IK,N), WMEAN(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST, RW, RWW, TEMP, TII, TJJ
      INTEGER           I, IERROR, IT, J, L, NREC
      LOGICAL           VARZ
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DCOPY, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
C     First test for errors in the input arguments
C
      IERROR = 1
      NREC = 1
      IF (MATRIX.NE.'V' .AND. MATRIX.NE.'v' .AND. MATRIX.NE.'R' .AND.
     *    MATRIX.NE.'r') THEN
         WRITE (P01REC,FMT=99999) MATRIX
      ELSE IF (K.LT.1) THEN
         WRITE (P01REC,FMT=99998) K
      ELSE IF (N.LT.2) THEN
         WRITE (P01REC,FMT=99997) N
      ELSE IF (IK.LT.K) THEN
         WRITE (P01REC,FMT=99996) IK, K
      ELSE IF (M.LT.1 .OR. M.GE.N) THEN
         WRITE (P01REC,FMT=99995) M, N
      ELSE
         IERROR = 0
         NREC = 0
C
C        Calculate mean of each series and the sample cross-covariance
C        matrix using a single pass update algorithm.
C
         DO 40 I = 1, K
            WMEAN(I) = W(I,1)
            DO 20 J = I, K
               R0(J,I) = ZERO
   20       CONTINUE
   40    CONTINUE
C
         DO 120 IT = 2, N
            RW = ONE/DBLE(IT)
            RWW = ONE - RW
            DO 80 I = 1, K
               DO 60 J = I, K
                  R0(J,I) = (WMEAN(I)-W(I,IT))*RWW*(WMEAN(J)-W(J,IT)) +
     *                      R0(J,I)
   60          CONTINUE
   80       CONTINUE
            DO 100 I = 1, K
               WMEAN(I) = (W(I,IT)-WMEAN(I))*RW + WMEAN(I)
  100       CONTINUE
C
  120    CONTINUE
C
C        Now calculate the cross-covariances at other lags.
C
         CONST = ONE/DBLE(N)
         DO 240 I = 1, K
            DO 220 J = 1, K
               CALL F06FBF(M,ZERO,R(I,J,1),IK*IK)
               DO 160 IT = 1, N - M
                  DO 140 L = 1, M
                     R(I,J,L) = R(I,J,L) + (W(I,IT)-WMEAN(I))*(W(J,IT+L)
     *                          -WMEAN(J))
  140             CONTINUE
  160          CONTINUE
               R(I,J,M) = R(I,J,M)*CONST
               DO 200 L = 1, M - 1
                  DO 180 IT = N - M + 1, N - L
                     R(I,J,L) = R(I,J,L) + (W(I,IT)-WMEAN(I))*(W(J,IT+L)
     *                          -WMEAN(J))
  180             CONTINUE
                  R(I,J,L) = R(I,J,L)*CONST
  200          CONTINUE
  220       CONTINUE
  240    CONTINUE
C
         IF (MATRIX.EQ.'V' .OR. MATRIX.EQ.'v') THEN
C
C           Copy the cross-covariances from the lower to the upper
C           triangle for lag zero.
C
            VARZ = .FALSE.
            DO 260 I = 1, K
               CALL DSCAL(K-I+1,CONST,R0(I,I),1)
               IF (I.LT.K) CALL DCOPY(K-I,R0(I+1,I),1,R0(I,I+1),IK)
               IF (R0(I,I).LE.ZERO) VARZ = .TRUE.
  260       CONTINUE
            IF (VARZ) THEN
               IERROR = 2
               NREC = 2
               WRITE (P01REC,FMT=99994)
            END IF
C
         ELSE
C
C           Convert cross-covariance matrices to cross-correlation
C           matrices. First convert diagonal elements at lag zero
C           to standard deviations checking that the varinances are
C           positive.
C
            VARZ = .FALSE.
            DO 320 I = 1, K
               R0(I,I) = R0(I,I)*CONST
               TEMP = R0(I,I)
               IF (TEMP.LE.ZERO) THEN
C
C                 Set elements involving this invalid variance to zero.
C
                  IF (I.GT.1) THEN
                     CALL F06FBF(I-1,ZERO,R0(I,1),IK)
                     DO 280 L = 1, M
                        CALL F06FBF(I-1,ZERO,R(I,1,L),IK)
                        CALL DCOPY(I-I,R(I,1,L),IK,R(1,I,L),1)
  280                CONTINUE
                  END IF
                  IF (I.LT.K) THEN
                     CALL F06FBF(K-I,ZERO,R0(I+1,I),1)
                     DO 300 L = 1, M
                        CALL F06FBF(K-I,ZERO,R(I+1,I,L),1)
                        CALL DCOPY(K-I,R(I+1,I,L),1,R(I,I+1,L),IK)
  300                CONTINUE
                  END IF
                  VARZ = .TRUE.
               ELSE
                  R0(I,I) = SQRT(TEMP)
               END IF
  320       CONTINUE
            IF (VARZ) THEN
               IERROR = 2
               NREC = 2
               WRITE (P01REC,FMT=99994)
C
C              There are variance errors so we must check the
C              elements R0(I,I).
C
               DO 360 J = 1, K - 1
                  TJJ = R0(J,J)
                  DO 340 I = J + 1, K
                     TII = R0(I,I)
                     IF (TII.GT.ZERO .AND. TJJ.GT.ZERO) R0(I,J) = R0(I,
     *                   J)*CONST/(TII*TJJ)
  340             CONTINUE
                  CALL DCOPY(K-J,R0(J+1,J),1,R0(J,J+1),IK)
  360          CONTINUE
               DO 420 L = 1, M
                  DO 400 J = 1, K
                     TJJ = R0(J,J)
                     DO 380 I = 1, K
                        TII = R0(I,I)
                        IF (TII.GT.ZERO .AND. TJJ.GT.ZERO) R(I,J,L)
     *                      = R(I,J,L)/(TII*TJJ)
  380                CONTINUE
  400             CONTINUE
  420          CONTINUE
            ELSE
C
C              There are no variance errors so dont need to check the
C              elements R0(I,I). Repeating above code but here there
C              are no IF controls, and this is the most likely route,
C              that is there will normally be no zero variances.
C
               DO 460 J = 1, K - 1
                  TJJ = R0(J,J)
                  DO 440 I = J + 1, K
                     R0(I,J) = R0(I,J)*CONST/(R0(I,I)*TJJ)
  440             CONTINUE
                  CALL DCOPY(K-J,R0(J+1,J),1,R0(J,J+1),IK)
  460          CONTINUE
               DO 520 L = 1, M
                  DO 500 J = 1, K
                     TJJ = R0(J,J)
                     DO 480 I = 1, K
                        R(I,J,L) = R(I,J,L)/(R0(I,I)*TJJ)
  480                CONTINUE
  500             CONTINUE
  520          CONTINUE
            END IF
         END IF
C
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT ('  ** On entry, MATRIX is an invalid character : MATRIX ',
     *       '= ',A1)
99998 FORMAT ('  ** On entry, K.lt.1 : K = ',I16)
99997 FORMAT ('  ** On entry, N.lt.2 : N = ',I16)
99996 FORMAT ('  ** On entry, IK.lt.K : IK = ',I16,' and K = ',I16)
99995 FORMAT ('  ** On entry, M.lt.1 or M.ge.N : M = ',I16,' and N = ',
     *       I16)
99994 FORMAT ('  ** On entry, at least one of the series is such that ',
     *       'all its elements are',/'     practically identical givin',
     *       'g zero (or near zero) variance')
      END
