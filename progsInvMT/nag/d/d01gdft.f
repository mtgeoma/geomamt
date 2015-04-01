      SUBROUTINE D01GDF(N,VECFUN,VECREG,NPTS,VK,NRAND,ITRANS,RES,ERR,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01GDF')
      INTEGER           IVSS, MAXDIM, MAXRUL
      PARAMETER         (IVSS=128,MAXDIM=20,MAXRUL=6)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERR, RES
      INTEGER           IFAIL, ITRANS, N, NPTS, NRAND
C     .. Array Arguments ..
      DOUBLE PRECISION  VK(N)
C     .. Subroutine Arguments ..
      EXTERNAL          VECFUN, VECREG
C     .. Local Scalars ..
      DOUBLE PRECISION  DPK, GRAD, PK, RES1, SUM, SUM2, XA, XX
      INTEGER           I, ICHECK, J, K, KUPPER, M, NDIM, NP, NREC, PTS
C     .. Local Arrays ..
      DOUBLE PRECISION  ALPHA(MAXDIM), C(IVSS), D(IVSS), FV(IVSS),
     *                  WT(IVSS), X(IVSS*MAXDIM)
      INTEGER           IVK(MAXDIM), KPTS(MAXRUL), KR(MAXRUL,MAXDIM)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MIN, MOD, NINT, DBLE, SQRT
C     .. Data statements ..
      DATA              KR(1,1), KR(2,1), KR(3,1), KR(4,1), KR(5,1),
     *                  KR(6,1)/1, 1, 1, 1, 1, 1/
      DATA              KR(1,2), KR(2,2), KR(3,2), KR(4,2), KR(5,2),
     *                  KR(6,2)/780, 1850, 3822, 6103, 15152, 30954/
      DATA              KR(1,3), KR(2,3), KR(3,3), KR(4,3), KR(5,3),
     *                  KR(6,3)/359, 1476, 544, 4104, 16592, 19394/
      DATA              KR(1,4), KR(2,4), KR(3,4), KR(4,4), KR(5,4),
     *                  KR(6,4)/766, 792, 1206, 6016, 9023, 15710/
      DATA              KR(1,5), KR(2,5), KR(3,5), KR(4,5), KR(5,5),
     *                  KR(6,5)/618, 840, 198, 6019, 12216, 2302/
      DATA              KR(1,6), KR(2,6), KR(3,6), KR(4,6), KR(5,6),
     *                  KR(6,6)/41, 2037, 2240, 4167, 4902, 9227/
      DATA              KR(1,7), KR(2,7), KR(3,7), KR(4,7), KR(5,7),
     *                  KR(6,7)/596, 229, 2304, 3851, 12506, 3420/
      DATA              KR(1,8), KR(2,8), KR(3,8), KR(4,8), KR(5,8),
     *                  KR(6,8)/86, 1578, 436, 4138, 7824, 3824/
      DATA              KR(1,9), KR(2,9), KR(3,9), KR(4,9), KR(5,9),
     *                  KR(6,9)/636, 526, 470, 259, 6093, 22300/
      DATA              KR(1,10), KR(2,10), KR(3,10), KR(4,10),
     *                  KR(5,10), KR(6,10)/287, 431, 1554, 1117, 12088,
     *                  5130/
      DATA              KR(1,11), KR(2,11), KR(3,11), KR(4,11),
     *                  KR(5,11), KR(6,11)/707, 1485, 480, 1188, 2399,
     *                  11222/
      DATA              KR(1,12), KR(2,12), KR(3,12), KR(4,12),
     *                  KR(5,12), KR(6,12)/707, 1450, 1004, 173, 8764,
     *                  17698/
      DATA              KR(1,13), KR(2,13), KR(3,13), KR(4,13),
     *                  KR(5,13), KR(6,13)/96, 1001, 684, 2919, 5491,
     *                  7057/
      DATA              KR(1,14), KR(2,14), KR(3,14), KR(4,14),
     *                  KR(5,14), KR(6,14)/49, 1001, 684, 235, 9274,
     *                  28739/
      DATA              KR(1,15), KR(2,15), KR(3,15), KR(4,15),
     *                  KR(5,15), KR(6,15)/373, 1001, 1447, 3043, 3054,
     *                  33207/
      DATA              KR(1,16), KR(2,16), KR(3,16), KR(4,16),
     *                  KR(5,16), KR(6,16)/613, 2, 857, 1249, 2648,
     *                  27717/
      DATA              KR(1,17), KR(2,17), KR(3,17), KR(4,17),
     *                  KR(5,17), KR(6,17)/373, 2, 2, 1249, 2648, 33207/
      DATA              KR(1,18), KR(2,18), KR(3,18), KR(4,18),
     *                  KR(5,18), KR(6,18)/2, 2, 2, 2, 2648, 1420/
      DATA              KR(1,19), KR(2,19), KR(3,19), KR(4,19),
     *                  KR(5,19), KR(6,19)/2, 2, 2, 2, 2, 1420/
      DATA              KR(1,20), KR(2,20), KR(3,20), KR(4,20),
     *                  KR(5,20), KR(6,20)/2, 2, 2, 2, 2, 2/
      DATA              KPTS(1), KPTS(2), KPTS(3), KPTS(4), KPTS(5),
     *                  KPTS(6)/2129, 5003, 10007, 20011, 40009, 80021/
C     .. Executable Statements ..
      NDIM = N
      PTS = NPTS
C     Validity check
      NREC = 1
      IF (NDIM.LT.1 .OR. NDIM.GT.MAXDIM) THEN
         WRITE (P01REC,FMT=99999) NDIM
         ICHECK = 1
         GO TO 260
      ELSE IF (NPTS.LT.1) THEN
         WRITE (P01REC,FMT=99998) NPTS
         ICHECK = 2
         GO TO 260
      ELSE IF (NRAND.LT.1) THEN
         WRITE (P01REC,FMT=99997) NRAND
         ICHECK = 3
         GO TO 260
      END IF
      NREC = 0
      ICHECK = 0
      IF (NPTS.LE.MAXRUL) THEN
C
C        Select Korobov vector
C
         NP = NPTS
         PTS = KPTS(NPTS)
         IVK(1) = 1
         IF (NDIM.EQ.1) GO TO 80
         IVK(2) = KR(NP,NDIM)
         IF (NDIM.EQ.2) GO TO 80
         XA = DBLE(IVK(2))
         DO 20 I = 3, NDIM
            IVK(I) = NINT(MOD(IVK(I-1)*XA,DBLE(PTS)))
   20    CONTINUE
         DO 40 I = 1, NDIM
            VK(I) = IVK(I)
   40    CONTINUE
      ELSE
         DO 60 I = 1, NDIM
            IVK(I) = NINT(VK(I))
   60    CONTINUE
      END IF
C
C     Begin integration
C
   80 DPK = 1.0D0/DBLE(PTS)
      SUM = 0.0D0
      SUM2 = 0.0D0
      DO 240 M = 1, NRAND
         RES1 = 0.0D0
C
C        Calculate random shift
C
         DO 100 K = 1, NDIM
            ALPHA(K) = G05CAF(ALPHA(K))
  100    CONTINUE
C
C        Calculate transformed integrand
C
         DO 220 I = 1, PTS, IVSS
            KUPPER = MIN(IVSS,PTS-I+1)
C
C           Initialise the weights
C
            DO 120 K = 1, KUPPER
               WT(K) = DPK
  120       CONTINUE
C
            DO 180 J = 1, NDIM
               CALL VECREG(NDIM,X,J,C,D,KUPPER)
               IF (ITRANS.NE.0) THEN
C
C                 No Periodizing transformation
C
                  DO 140 K = 1, KUPPER
                     PK = DBLE(I+K-1)
                     XX = ALPHA(J) + IVK(J)*PK*DPK
                     XX = XX - INT(XX)
                     GRAD = D(K) - C(K)
                     WT(K) = WT(K)*GRAD
                     X(K+(J-1)*KUPPER) = C(K) + GRAD*XX
  140             CONTINUE
               ELSE IF (ITRANS.EQ.0) THEN
C
C                 Periodizing transformation
C
                  DO 160 K = 1, KUPPER
                     PK = DBLE(I+K-1)
                     XX = ALPHA(J) + IVK(J)*PK*DPK
                     XX = XX - INT(XX)
                     GRAD = D(K) - C(K)
                     WT(K) = WT(K)*GRAD
                     WT(K) = WT(K)*6.0D0*XX*(1.0D0-XX)
                     XX = XX*XX*(3.0D0-2.0D0*XX)
                     X(K+(J-1)*KUPPER) = C(K) + GRAD*XX
  160             CONTINUE
               END IF
  180       CONTINUE
            CALL VECFUN(NDIM,X,FV,KUPPER)
            DO 200 K = 1, KUPPER
               RES1 = RES1 + FV(K)*WT(K)
  200       CONTINUE
  220    CONTINUE
         SUM = SUM + RES1
         SUM2 = SUM2 + RES1*RES1
  240 CONTINUE
C
C     Calculate the mean
C
      RES = SUM/DBLE(NRAND)
      IF (NRAND.EQ.1) THEN
         ERR = 0.0D0
      ELSE
C
C        Standard error
C
         ERR = SQRT(ABS((SUM2-DBLE(NRAND)*RES*RES)/DBLE(NRAND*(NRAND-1))
     *         ))
      END IF
  260 IFAIL = P01ABF(IFAIL,ICHECK,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, NDIM must satisfy 1 .LE. NDIM .LE. 20: ND',
     *       'IM = ',I12)
99998 FORMAT (' ** On entry, NPTS must be at least 1: NPTS = ',I12)
99997 FORMAT (' ** On entry, NRAND must be at least 1: NRAND = ',I12)
      END
