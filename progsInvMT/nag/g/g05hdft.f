      SUBROUTINE G05HDF(MODE,K,IP,IQ,MEAN,PAR,LPAR,QQ,IK,N,W,REF,LREF,
     *                  IWORK,LIWORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1682 (JUN 1995).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05HDF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, IP, IQ, K, LIWORK, LPAR, LREF, N
      CHARACTER         MEAN, MODE
C     .. Array Arguments ..
      DOUBLE PRECISION  PAR(LPAR), QQ(IK,K), REF(LREF), W(IK,N)
      INTEGER           IWORK(LIWORK)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGEST, EPS, LIMIT, SUM2
      INTEGER           I, IERR, IERROR, IFAILA, IPK, IQK, IR, J, K3,
     *                  KP, KQ, KW, L, L3, LC, LG, LM, LR2, LR3, LW1,
     *                  LW2, LW3, LW4, LW5, LW6, LW7, LW8, LW9, M, MK,
     *                  NPAR, NREC, PK2, T
      LOGICAL           STAT
      CHARACTER         UMEAN, UMODE
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04UDU, F03ABF, DAXPY, DCOPY, F06FBF, DGEMV,
     *                  G05HDU, G05HDV, G05HDW, G05HDX, G05HDZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, SQRT
C     .. Executable Statements ..
C
C     ******************************************************************
C
C     This routine generates a realization of the recent history of a
C     multivariate time series which is stored in the array REF. The
C     array REF is also used to store the details of the spectral
C     decompositions of GAMMA(0) and SIGMA. This allows a new recent
C     history to be generated using MODE = 'Restart', without having to
C     recompute these spectral decompositions.
C     The routine may also be called with MODE = 'Continue' when a
C     further sequence is required from the point at which the last
C     sequence ended. This avoids the recomputing of both the spectral
C     decompositions and the recent history.
C
C     The computational procedure is described in a paper by Shea, B.L.
C        'A note on the generation of independent realizations of a
C         vector autoregressive moving average process'
C         Journal of Time Series Analysis, Vol 9, pages 403 - 410, 1988.
C     This paper however advocates the use of a Cholesky factorisation.
C     Although the spectral decomposition is roughly 10-20 percent
C     slower it does give more stable results on a range of computers.
C     Thus instead of using the current NAG routines (G05EAF and G05EZF)
C     to generate the multivarite Normals (using Cholesky), this routine
C     uses two auxillary routines adpated from these G05 routines to use
C     the spectral decomposition.
C
C     ******************************************************************
C
C     First test for errors in the input arguments.
C
      IERROR = 0
      NREC = 0
      UMODE = MODE
      UMEAN = MEAN
      CALL E04UDU(UMODE)
      CALL E04UDU(UMEAN)
      NPAR = (IP+IQ)*K*K
      IF (UMEAN.EQ.'M') NPAR = NPAR + K
      IERROR = 1
      NREC = 1
      IF (UMODE.NE.'S' .AND. UMODE.NE.'R' .AND. UMODE.NE.'C') THEN
         WRITE (P01REC,FMT=99999) MODE
      ELSE IF (UMEAN.NE.'M' .AND. UMEAN.NE.'Z') THEN
         WRITE (P01REC,FMT=99998) MEAN
      ELSE IF (K.LT.1) THEN
         WRITE (P01REC,FMT=99997) K
      ELSE IF (IP.LT.0) THEN
         WRITE (P01REC,FMT=99996) IP
      ELSE IF (IQ.LT.0) THEN
         WRITE (P01REC,FMT=99995) IQ
      ELSE IF (LPAR.LT.NPAR) THEN
         NREC = 2
         WRITE (P01REC,FMT=99994) LPAR, NPAR
      ELSE IF (N.LT.1) THEN
         WRITE (P01REC,FMT=99993) N
      ELSE IF (IK.LT.K) THEN
         WRITE (P01REC,FMT=99992) IK, K
      ELSE
         IERROR = 0
         NREC = 0
      END IF
      IF (IERROR.GT.0) GO TO 420
C
      IF (K.GE.6) THEN
         M = K*MAX(IP,IQ)
      ELSE
         M = K*(IP+IQ)
      END IF
      LG = (M+1)*(M+2)
      LC = (K+1)*(K+2)
      LR2 = M + LC
      LR3 = LR2
      IF (IP.NE.0) LR3 = LR3 + LG
C
C     Test whether the REF array is big enough.
C
      L = K*(K+1)/2
      IF (IP.GE.2) L = L + (IP-1)*K*K
      J = M*M + L*(L+3) + K*K*(IQ+1)
      I = K*MAX(IP,IQ)*(K*MAX(IP,IQ)+2)
      IF (K.GE.6) THEN
         J = 4*M*M
      ELSE
         J = MAX(I,J)
      END IF
      J = LR3 + J
C
      IF (LREF.LT.J) THEN
         IERROR = 1
         NREC = 2
         WRITE (P01REC,FMT=99991) LREF, J
         GO TO 420
      END IF
C
C     Test whether the integer workspace array is big enough.
C
      J = K*MAX(IP,IQ)
      IF (LIWORK.LT.J) THEN
         IERROR = 1
         NREC = 2
         WRITE (P01REC,FMT=99990) LIWORK, J
         GO TO 420
      END IF
C
C     Only for first call.
C
      IF (UMODE.EQ.'S') THEN
C
C        Test whether QQ is positive-definite.
C
C        First set the upper triangle of QQ to the lower triangle
C        of QQ.
C
         DO 20 I = 1, K - 1
            CALL DCOPY(K-I,QQ(I+1,I),1,QQ(I,I+1),IK)
   20    CONTINUE
C
         IFAILA = 1
         CALL F03ABF(QQ,IK,K,SUM2,REF,IFAILA)
C
C        Reset the strict lower triangle of QQ (which F03ABF has
C        corrupted) to its former value.
C
         DO 40 I = 1, K - 1
            CALL DCOPY(K-I,QQ(I,I+1),IK,QQ(I+1,I),1)
   40    CONTINUE
C
         IF (IFAILA.GT.0) THEN
            IERROR = 2
            NREC = 1
            WRITE (P01REC,FMT=99989)
            GO TO 420
         END IF
C
C        Test for stationarity of the AR operator.
C
         IF (IP.GT.0) THEN
            KP = K*IP
            LW1 = LR3 + 1
            LW2 = LW1 + KP*KP
            LW3 = LW2 + KP
            LIMIT = 1.0D0
            CALL G05HDZ(IP,K,PAR,REF(LW1),REF(LW2),REF(LW3),IWORK,KP,
     *                  LIMIT,BIGEST,STAT,IERR)
            IF (IERR.EQ.1) THEN
               IERROR = 3
               NREC = 2
               WRITE (P01REC,FMT=99988)
               GO TO 420
            ELSE IF ( .NOT. STAT) THEN
               IERROR = 2
               NREC = 1
               WRITE (P01REC,FMT=99987)
               GO TO 420
            END IF
         END IF
C
C        Test for invertibility of the MA operator.
C
         IF (IQ.GT.0) THEN
            LW4 = IP*K*K + 1
            KQ = K*IQ
            LW1 = LR3 + 1
            LW2 = LW1 + KQ*KQ
            LW3 = LW2 + KQ
            LIMIT = 1.0D0 + SQRT(X02AJF())
            CALL G05HDZ(IQ,K,PAR(LW4),REF(LW1),REF(LW2),REF(LW3),IWORK,
     *                  KQ,LIMIT,BIGEST,STAT,IERR)
            IF (IERR.EQ.1) THEN
               IERROR = 3
               NREC = 2
               WRITE (P01REC,FMT=99988)
               GO TO 420
            ELSE IF ( .NOT. STAT) THEN
               IERROR = 2
               NREC = 1
               WRITE (P01REC,FMT=99986)
               GO TO 420
            END IF
         END IF
C
C        Set up first part of reference vector REF (for SIGMA matrix).
C
         CALL F06FBF(K,0.0D0,REF(LR3+1),1)
         EPS = 2.0D0*X02AJF()
         IFAILA = 1
         CALL G05HDV(REF(LR3+1),K,QQ,IK,EPS,REF(M+1),LC,IFAILA)
         IF (IFAILA.GT.0) THEN
            IERROR = 4
            NREC = 2
            WRITE (P01REC,FMT=99985)
            GO TO 420
         END IF
      END IF
C
C     First deal with pure MA models.
C
      IF (IP.EQ.0) THEN
         IF (UMODE.NE.'C') THEN
            LW2 = 1
            DO 60 I = 1, IQ
               IFAILA = 1
               CALL G05HDU(REF(LW2),K,REF(M+1),LC,IFAILA)
               IF (IFAILA.GT.0) THEN
                  IERROR = 6
                  NREC = 2
                  WRITE (P01REC,FMT=99984)
                  GO TO 420
               END IF
               LW2 = LW2 + K
   60       CONTINUE
         END IF
         GO TO 80
      END IF
C
      IF (UMODE.EQ.'S') THEN
C
C        Now deal with ARMA and pure AR models.
C
         IF ((K.GE.6) .OR. (IQ.EQ.0 .AND. K.GE.5)) THEN
C
C           Use the doubling algorithm to compute GAMMA(0).
C
            LW3 = LR3 + 2*M*M + 1
            CALL G05HDW(K,IP,IQ,PAR,NPAR,QQ,IK,REF(LW3),M,REF(LR3+1),
     *                  IERROR)
C
         ELSE
C
C           Use the direct method of computing GAMMA(0) with Barone's
C           state space representation.
C
            KW = K*(K+1)/2
            IF (IP.GE.2) KW = KW + (IP-1)*K*K
            LW5 = LR3 + M*M + 1
            LW6 = LW5 + KW
            LW7 = LW6 + KW
            LW8 = LW7 + (IQ+1)*K*K
            LW9 = LW8 + KW*KW
            CALL G05HDX(K,IP,IQ,PAR,NPAR,M,QQ,IK,REF(LR3+1),KW,REF(LW5),
     *                  REF(LW6),REF(LW7),REF(LW8),REF(LW9),IERROR)
         END IF
C
         IF (IERROR.GT.0) THEN
            IERROR = 5
            NREC = 2
            WRITE (P01REC,FMT=99983)
            GO TO 420
         END IF
C
C        Set up rest of reference vector.
C
         CALL F06FBF(M,0.0D0,REF,1)
         EPS = 2.0D0*X02AJF()
         IFAILA = 1
         CALL G05HDV(REF,M,REF(LR3+1),M,EPS,REF(LR2+1),LG,IFAILA)
         IF (IFAILA.GT.0) THEN
            IERROR = 4
            NREC = 2
            WRITE (P01REC,FMT=99981)
            GO TO 420
         END IF
      END IF
C
      IF (UMODE.NE.'C') THEN
         IFAILA = 1
         CALL G05HDU(REF,M,REF(LR2+1),LG,IFAILA)
         IF (IFAILA.GT.0) THEN
            IERROR = 6
            NREC = 2
            WRITE (P01REC,FMT=99984)
            GO TO 420
         END IF
      END IF
C
C     ***************************************************************
C
C     This section generates the realization of the multivariate time
C     series. A realization of the recent history of the series is
C     stored in the array REF and updated on exit. A user would
C     normally summon this section directly using MODE = 'Continue' if
C     they have already called G05HDF.
C
   80 CONTINUE
      K3 = K*K
      PK2 = IP*K3
      IF ((K.LE.5) .OR. (IP.EQ.0) .OR. (IQ.EQ.0)) THEN
C
         KP = K*IP
         DO 180 T = 1, N
C
C           Generate a(t).
C
            IFAILA = 1
            CALL G05HDU(W(1,T),K,REF(M+1),LC,IFAILA)
            CALL DCOPY(K,W(1,T),1,REF(LR3+1),1)
            IF (IFAILA.GT.0) THEN
               IERROR = 6
               NREC = 1
               WRITE (P01REC,FMT=99982)
               GO TO 420
            END IF
C
C           Calculate z(t) from the ARMA equation.
C
            DO 100 L = 1, IP
               CALL DGEMV('T',K,K,1.0D0,PAR((L-1)*K3+1),K,REF((L-1)
     *                     *K+1),1,1.0D0,W(1,T),1)
  100       CONTINUE
            DO 120 L = 1, IQ
               CALL DGEMV('T',K,K,-1.0D0,PAR(PK2+(L-1)*K3+1),K,
     *                     REF(KP+(L-1)*K+1),1,1.0D0,W(1,T),1)
  120       CONTINUE
C
C           Update REF array.
C
            DO 140 L = IP, 2, -1
               CALL DCOPY(K,REF((L-2)*K+1),1,REF((L-1)*K+1),1)
  140       CONTINUE
C
            DO 160 L = IQ, 2, -1
               CALL DCOPY(K,REF(KP+(L-2)*K+1),1,REF(KP+(L-1)*K+1),1)
  160       CONTINUE
C
            IF (IP.GT.0) CALL DCOPY(K,W(1,T),1,REF,1)
            IF (IQ.GT.0) CALL DCOPY(K,REF(LR3+1),1,REF(KP+1),1)
C
  180    CONTINUE
C
      ELSE
C
C        Generate a(0).
C
         IFAILA = 1
         CALL G05HDU(W(1,1),K,REF(M+1),LC,IFAILA)
         IF (IFAILA.GT.0) THEN
            IERROR = 6
            NREC = 1
            WRITE (P01REC,FMT=99982)
            GO TO 420
         END IF
         CALL DCOPY(K,W(1,1),1,REF(LR3+1),1)
C
C        Put z(0) = h' * ALPHA(0)  +  a(0)
C
         CALL DAXPY(K,1.0D0,REF,1,W(1,1),1)
C
C        Set REF(LR3+1) = R
C
         IR = MAX(IP,IQ)
         LM = K + IR*K*K
         DO 220 L = 1, IR
            DO 200 I = 1, K
               CALL F06FBF(K,0.0D0,REF(LR3+K+(L-1)*K3+(I-1)*K+1),1)
  200       CONTINUE
  220    CONTINUE
C
         DO 260 L = 1, IP
            DO 240 I = 1, K
               CALL DCOPY(K,PAR((L-1)*K3+(I-1)*K+1),1,REF(LR3+K+(L-1)
     *                     *K3+(I-1)*K+1),1)
  240       CONTINUE
  260    CONTINUE
C
         DO 300 L = 1, IQ
            DO 280 I = 1, K
               CALL DAXPY(K,-1.0D0,PAR(PK2+(L-1)*K3+(I-1)*K+1),1,
     *                     REF(LR3+K+(L-1)*K3+(I-1)*K+1),1)
  280       CONTINUE
  300    CONTINUE
C
         IF (IP.LE.IQ) THEN
C
            DO 340 T = 1, N
C
C              Set ALPHA(t) = T * ALPHA(t-1)  +  R * a(t-1)
C
               MK = IR*K
               CALL DGEMV('T',K,MK,1.0D0,REF(LR3+K+1),K,REF(LR3+1),1,
     *                     0.0D0,REF(LR3+LM+1),1)
C
               IPK = IP*K
               CALL DGEMV('T',K,IPK,1.0D0,PAR,K,REF,1,1.0D0,
     *                     REF(LR3+LM+1),1)
C
               DO 320 L = 1, IR - 1
                  CALL DAXPY(K,1.0D0,REF(L*K+1),1,REF(LR3+LM+(L-1)*K+1)
     *                        ,1)
  320          CONTINUE
C
C              Update ALPHA(t)
C
               CALL DCOPY(MK,REF(LR3+LM+1),1,REF,1)
C
               IF (T.LT.N) THEN
C
C                 Generate a(t)
C
                  IFAILA = 1
                  CALL G05HDU(W(1,T+1),K,REF(M+1),LC,IFAILA)
                  IF (IFAILA.GT.0) THEN
                     IERROR = 6
                     NREC = 1
                     WRITE (P01REC,FMT=99982)
                     GO TO 420
                  END IF
                  CALL DCOPY(K,W(1,T+1),1,REF(LR3+1),1)
C
C                 Put z(t) = h' * ALPHA(t)  +  a(t)
C
                  CALL DAXPY(K,1.0D0,REF,1,W(1,T+1),1)
C
               END IF
C
  340       CONTINUE
C
         ELSE
C
            DO 380 T = 1, N
C
C              Set ALPHA(t) = A * ALPHA(t-1)  +  R * z(t-1)
C
               MK = IR*K
               CALL DGEMV('T',K,MK,1.0D0,REF(LR3+K+1),K,W(1,T),1,0.0D0,
     *                     REF(LR3+LM+1),1)
C
               IQK = IQ*K
               CALL DGEMV('T',K,IQK,1.0D0,PAR(PK2+1),K,REF,1,1.0D0,
     *                     REF(LR3+LM+1),1)
C
               DO 360 L = 1, IR - 1
                  CALL DAXPY(K,1.0D0,REF(L*K+1),1,REF(LR3+LM+(L-1)*K+1)
     *                        ,1)
  360          CONTINUE
C
C              Update ALPHA(t)
C
               CALL DCOPY(MK,REF(LR3+LM+1),1,REF,1)
C
               IF (T.LT.N) THEN
C
C                 Generate a(t)
C
                  IFAILA = 1
                  CALL G05HDU(W(1,T+1),K,REF(M+1),LC,IFAILA)
                  IF (IFAILA.GT.0) THEN
                     IERROR = 6
                     NREC = 1
                     WRITE (P01REC,FMT=99982)
                     GO TO 420
                  END IF
C
C                 Put z(t) = h' * ALPHA(t)  +  a(t)
C
                  CALL DAXPY(K,1.0D0,REF,1,W(1,T+1),1)
C
               END IF
C
  380       CONTINUE
C
         END IF
C
C        Add the mean on to z(t)
C
         IF (UMEAN.EQ.'M') THEN
            L3 = (IP+IQ)*K3
            DO 400 T = 1, N
               CALL DAXPY(K,1.0D0,PAR(L3+1),1,W(1,T),1)
  400       CONTINUE
         END IF
C
      END IF
C
  420 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, MODE is an invalid character : MODE = ',A1)
99998 FORMAT (' ** On entry, MEAN is an invalid character : MEAN = ',A1)
99997 FORMAT (' ** On entry, K.lt.1 : K = ',I16)
99996 FORMAT (' ** On entry, IP.lt.0 : IP = ',I16)
99995 FORMAT (' ** On entry, IQ.lt.0 : IQ = ',I16)
99994 FORMAT (' ** On entry, LPAR .lt. the number of parameters : LPAR',
     *       ' = ',I16,/'    and the number of parameters = ',I16)
99993 FORMAT (' ** On entry, N.lt.1 : N = ',I16)
99992 FORMAT (' ** On entry, IK.lt.K : IK = ',I16,'  and K = ',I16)
99991 FORMAT (' ** On entry, LREF is too small : LREF = ',I16,'  but m',
     *       'ust be at',/'    least ',I16)
99990 FORMAT (' ** On entry, LIWORK is too small : LIWORK = ',I16,'  b',
     *       'ut must be at',/'    least ',I16)
99989 FORMAT (' ** On entry, the covariance matrix QQ is not positive-',
     *       'definite')
99988 FORMAT (' ** An excessive number of iterations were required by ',
     *       'F02AFF to evaluate ',/'   the eigenvalues of the matrice',
     *       's used to test for stationarity or invertibility.')
99987 FORMAT (' ** On entry, the AR parameter matrices are outside the',
     *       ' stationarity region')
99986 FORMAT (' ** On entry, the MA parameter matrices are outside the',
     *       ' invertibility region')
99985 FORMAT (' ** An excessive number of iterations were required by ',
     *       'F02ABF to evaluate ',/'   the eigenvalues of the covaria',
     *       'nce matrix.')
99984 FORMAT (' ** On entry, MODE was set to R or r but it has not bee',
     *       'n possible to',/'    generate a realisation of the recen',
     *       't history of the series')
99983 FORMAT (' ** The reference vector cannot be computed because the',
     *       ' AR parameters are too',/'     close to the boundary of ',
     *       'the stationarity region')
99982 FORMAT (' ** On entry the reference vector REF has been corrupte',
     *       'd between calls to G05HDF')
99981 FORMAT (' ** An excessive number of iterations were required by ',
     *       'F02ABF to evaluate ',/'   the eigenvalues to be stored i',
     *       'n the reference vector.')
      END
