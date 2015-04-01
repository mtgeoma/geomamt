      SUBROUTINE G07EAW(X,N,K1,K2,LOWER,CHECKP,A,PROB,EPS,IWRK,WRK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C
C     Computes either the Hodges-Lehmann location estimator:
C        median of (X(I) + X(J))/2   for 1.le.I.le.J.le.N
C     or one of the confidence interval limits.
C
C     Arguments  X   Real array of observations
C                 * Values of X must be in nondecreasing order *
C
C                N   INTEGER number of observations  (INPUT)
C
C                K1 and K2  -  define which ordered average to find
C
C                LOWER  LOGICAL for lower or upper tail
C
C                CHECKP  LOGICAL - TRUE for confidence limits which
C                        require probability estimates, FALSE for the
C                        central estimate
C
C                A  REAL the final estimate
C
C                PROB  REAL the tail probability for one of the
C                      confidence limits
C
C                EPS  a small multiple of machine precision
C
C                IWRK         INTEGER array of length N for workspace
C
C                IWRK(N+1)    INTEGER array of length N for workspace
C
C                IWRK(2*N+1)  INTEGER array of length N for workspace
C
C     External routine
C        G05CAF  Function providing uniform random variables
C                in the interval (0,1).
C
C     Notes
C        G07EAW has an expected time complexity on the
C        order of N * LOG(N)
C
C     Based on the Function HLQEST.
C     J F Monahan, April 1982, Dept of Stat, NCSU, Raleigh, NC 27650
C     Final version  June 1983
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, EPS, PROB
      INTEGER           K1, K2, N
      LOGICAL           CHECKP, LOWER
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(3*N), X(N)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AM, AMN, AMX, W, WN, XMED
      INTEGER           I, IF2, IPIQ, IQ, J, K, L, M, MDLL, MDLROW,
     *                  MDLU, N1, N2, SM, SQ
      LOGICAL           COND, DONE
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF, X02AJF
      EXTERNAL          G05CAF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          G07EAU, G07EAV, G08AGF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, MIN
C     .. Executable Statements ..
      M = N*(N+1)/2
      N2 = 2*N
C
C     Take care of special case: N=2
C
      IF (N.EQ.2) THEN
         IF (K1.EQ.3) A = X(2)
         IF (K1.EQ.2) A = (X(1)+X(2))/2.D0
         IF (K1.EQ.1) A = X(1)
         IF (K1.EQ.0) A = X(1)
      ELSE
C
C        Find the total number of pairs (M) and the
C        median(s) (K1,K2) needed
C
C        Initialize left and right bounds
C
         DO 20 I = 1, N
            IWRK(I) = I
            IWRK(N+I) = N
   20    CONTINUE
C        SM = Number in set S at step M
         SM = M
C        L = Number of pairs less than those in set S at step M
         L = -1
C
C        Use the median of X(I)'s to partition on the first step
C
         MDLL = K1*2/(N+1)
         MDLU = K2*2/(N+1)
         IF (MDLL.LT.1) MDLL = 1
         IF (MDLU.LT.1) MDLU = 1
         AM = X(MDLL) + X(MDLU)
         AMN = X(1) + X(1)
         AMX = X(N) + X(N)
         DONE = .FALSE.
   40    CONTINUE
         COND = ( .NOT. DONE) .AND. (ABS(AMN-AMX).GT.(2.0D0*X02AJF()))
     *           .AND. (SM.NE.2 .OR. (SM.EQ.2 .AND. K1.EQ.K2))
C
C        Unless finished start with partition step
C
         IF (COND) THEN
C
C           Go and partition.
C
            CALL G07EAV(AM,N,X,IWRK,SQ)
C
C           ***  Finished partition --- Start branching  ***
C
            IF (SQ.EQ.L) THEN
C
C              If consecutive partitions are the same we probably
C              have ties
C
               CALL G07EAU(N,X,IWRK,AM,AMN,AMX)
            ELSE
C
C              Are we nearly done, with the values we want on the
C              border?
C
               IF (SQ.EQ.K2-1 .OR. SQ.EQ.K1) THEN
C
C
C                 Find:   Max of those .lt. AM
C                         Min of those .ge. AM
C
                  AMN = X(N) + X(N)
                  AMX = X(1) + X(1)
                  DO 60 I = 1, N
                     IQ = IWRK(N2+I)
                     IPIQ = I + IQ
                     IF (IQ.GT.0) AMX = MAX(AMX,X(I)+X(IPIQ-1))
                     IPIQ = I + IQ
                     IF (IQ.LT.N-I+1) AMN = MIN(AMN,X(I)+X(IPIQ))
   60             CONTINUE
                  AM = (AMN+AMX)/2.0D0
C
C                 We are done, but which situation are we in?
C
                  IF (K1.GE.K2) THEN
                     IF (SQ.EQ.K1) AM = AMX
                     IF (SQ.EQ.K1-1) AM = AMN
                  END IF
                  DONE = .TRUE.
               ELSE
C
C                 The set S is split, which piece do we keep?
C
                  IF (SQ.GT.K1) THEN
C
C                    New S = (Old S) .INTERSECT. (Those .lt. AM)
C                    Cut off top
C
                     DO 80 I = 1, N
C
C                       Reset right bounds for each row
C
                        IWRK(N+I) = I + IWRK(N2+I) - 1
   80                CONTINUE
                  ELSE
C
C                    Cut off bottom
C                    New S = (Old S) .INTERSECT. (Those .ge. AM)
C
                     DO 100 I = 1, N
C
C                       Reset left bounds for each row
C
                        IWRK(I) = I + IWRK(N2+I)
  100                CONTINUE
                  END IF
C
C                 Count  SM = Number of pairs still in new set S
C                  L  = Number of pairs less than those in new set S
C
                  L = 0
                  SM = 0
                  DO 120 I = 1, N
                     L = L + IWRK(I) - I
                     SM = SM + IWRK(N+I) - IWRK(I) + 1
  120             CONTINUE
C
C                 ***  Normal restart jump  ***
C
C
                  IF (SM.LE.2) THEN
C
C                    IF SM.le.2 then can only get to 2 left if K1.ne.K2
C                    Go get their average
C
                     CALL G07EAU(N,X,IWRK,AM,AMN,AMX)
C
                  ELSE
C
C                    ***  Restart here unless worried about ties  ***
C
C                    Use random row median as partition element
C
                     K = INT(DBLE(SM)*G05CAF(0.0D0))
C
C                    K is a random integer from 0 to SM-1
C
                     DO 140 I = 1, N
                        J = I
                        IF (K.LE.IWRK(N+I)-IWRK(I)) THEN
                           GO TO 160
                        ELSE
                           K = K - IWRK(N+I) + IWRK(I) - 1
                        END IF
  140                CONTINUE
C
C                    J is a random row -- now get its median
C
  160                MDLROW = (IWRK(J)+IWRK(N+J))/2
                     AM = X(J) + X(MDLROW)
                  END IF
               END IF
            END IF
C
C           Return to top of WHILE COND DO loop
C
            GO TO 40
         END IF
         A = AM/2.0D0
      END IF
C
C     For the confidence limits use the Wilcoxon statistic adjusted
C     for ties to check that the tail probability is correct.
C
      IF (CHECKP) THEN
         IF (LOWER) THEN
            XMED = A - EPS
            CALL G08AGF(N,X,XMED,'UPPER','Y',W,WN,PROB,N1,WRK,IF2)
         ELSE
            XMED = A + EPS
            CALL G08AGF(N,X,XMED,'LOWER','Y',W,WN,PROB,N1,WRK,IF2)
         END IF
      END IF
      RETURN
      END
