      SUBROUTINE G07EBW(N,X,M,Y,K1,K2,A,IWRK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Purpose    Find value of A = Y(J) - X(I) for 1.le.I.le.N and
C                1.le.J.le.M corresponding to given K1 and K2.
C                Adapted from routine to compute the Hodges-Lehmann
C                location estimator for one sample.
C     Based on the Function HLQEST.
C     J F Monahan, April 1982, Dept of Stat, NCSU, Raleigh, NC 27650
C     Final version  June 1983
C
C     Arguments  X   Real array of N observations  (INPUT)
C                 * Values of X must be in nondecreasing order *
C                Y   Real array of M observations  (INPUT)
C                 * Values of Y must be in nondecreasing order *
C
C                K1 and K2
C
C                A           The result
C
C                PROB        Estimate of tail probability if CHECKP true
C
C                IWRK        INTEGER array of length 3*N for workspace
C
C     External routine
C        G05CAF  Function providing uniform random variables
C                in the interval (0,1).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A
      INTEGER           K1, K2, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(M)
      INTEGER           IWRK(3*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AM, AMN, AMX
      INTEGER           I, IPIQ, IQ, IX, IY, J, K, L, MDLROW, N2, NM,
     *                  SM, SQ
      LOGICAL           COND, DONE
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      EXTERNAL          G05CAF
C     .. External Subroutines ..
      EXTERNAL          G07EBU, G07EBV
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN
C     .. Executable Statements ..
C
C     Take care of special case: N=2 and M=2
C
C     Find the total number of pairs (NM) and the
C
      NM = N*M
      N2 = 2*N
C
C     Initialize left and right bounds
C
      DO 20 I = 1, N
         IWRK(I) = 1
         IWRK(N+I) = M
   20 CONTINUE
C
C     SM = Number in set S at step M
C
      SM = NM
C
C     L = Number of pairs less than those in set S at step M
C
      L = -1
C
C     Use rough estimates based on K1 and K2 to partition on
C     the first step
C
      IY = MAX(1,K2/N)
      IX = MAX(1,K1/M)
      AM = Y(IY) - X(N-IX+1)
      AMN = Y(1) - X(N)
      AMX = Y(M) - X(1)
      DONE = .FALSE.
   40 CONTINUE
      COND = ( .NOT. DONE) .AND. (AMN.NE.AMX .AND. SM.NE.2)
C
C     Unless finished start with partition step
C
      IF (COND) THEN
C
C        Go and partition.
C
         CALL G07EBV(AM,N,X,M,Y,IWRK,SQ)
C
C        ***  Finished partition --- Start branching  ***
C
         IF (SQ.EQ.L) THEN
C
C           If consecutive partitions are the same we probably have ties
C
            CALL G07EBU(N,X,M,Y,IWRK,AM,AMN,AMX)
         ELSE
C
C           Are we nearly done, with the values we want on the border?
C           IF (We need  max of those .Lt. AM -or- min of those .ge. AM)
C
            IF (SQ.EQ.K2-1 .OR. SQ.EQ.K1) THEN
C
C              Find:   Max of those .lt. AM
C                      Min of those .ge. AM
C
               AMN = Y(M) - X(1)
               AMX = Y(1) - X(N)
               DO 60 I = 1, N
                  IQ = IWRK(N2+I)
                  IPIQ = IQ + 1
                  IF (IQ.GT.0) AMX = MAX(AMX,Y(IPIQ-1)-X(N-I+1))
                  IPIQ = IQ + 1
                  IF (IQ.LT.M) AMN = MIN(AMN,Y(IPIQ)-X(N-I+1))
   60          CONTINUE
               AM = (AMN+AMX)/2.0D0
C
C              We are done, but which situation are we in?
C
               IF (K1.GE.K2) THEN
                  IF (SQ.EQ.K1) AM = AMX
                  IF (SQ.EQ.K1-1) AM = AMN
               END IF
               DONE = .TRUE.
            ELSE
C
C              The set S is split, which piece do we keep?
C
               IF (SQ.GT.K1) THEN
C
C                 New S = (Old S) .INTERSECT. (Those .lt. AM)
C                 Cut off top
C
                  DO 80 I = 1, N
C
C                    Reset right bounds for each row
C
                     IWRK(N+I) = IWRK(N2+I)
   80             CONTINUE
               ELSE
C
C                 Cut off bottom
C                 New S = (Old S) .INTERSECT. (Those .ge. AM)
C
                  DO 100 I = 1, N
C
C                    Reset left bounds for each row
C
                     IWRK(I) = IWRK(N2+I) + 1
  100             CONTINUE
               END IF
C
C              Count  SM = Number of pairs still in new set S
C                     L  = Number of pairs less than those in new set S
C
               L = 0
               SM = 0
               DO 120 I = 1, N
                  L = L + IWRK(I) - 1
                  SM = SM + IWRK(N+I) - IWRK(I) + 1
  120          CONTINUE
C
C              ***  Normal restart jump  ***
C
               IF (SM.LE.2) THEN
C
C                 IF SM.le.2 then can only get to 2 left if K1.ne.K2
C                 Go get their average
C
                  CALL G07EBU(N,X,M,Y,IWRK,AM,AMN,AMX)
C
               ELSE
C
C                 ***  Restart here unless worried about ties  ***
C
C                 Use random row median as partition element
C
                  K = INT(DBLE(SM)*G05CAF(0.0D0))
C
C                 K is a random integer from 0 to SM-1
C
                  DO 140 I = 1, N
                     J = I
                     IF (K.LE.IWRK(N+I)-IWRK(I)) THEN
                        GO TO 160
                     ELSE
                        K = K - IWRK(N+I) + IWRK(I) - 1
                     END IF
  140             CONTINUE
C
C                 J is a random row -- now get its median
C
  160             MDLROW = (IWRK(J)+IWRK(N+J))/2
                  AM = Y(MDLROW) - X(N-J+1)
               END IF
            END IF
         END IF
C
C        Return to top of WHILE COND DO loop
C
         GO TO 40
      END IF
      A = AM
      RETURN
      END
