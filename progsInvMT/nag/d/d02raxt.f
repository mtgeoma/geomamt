      SUBROUTINE D02RAX(K,P,Q,N,M,A,X,NMAX,Y,S,MTNMAX,ALF,C,IERROR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 10 REVISED. IER-377 (JUN 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     *******************************************************
C     *                                                     *
C     *   THIS IS A TWO POINT BOUNDARY VALUE DEFERRED       *
C     *   CORRECTION GENERATOR FOR SYSTEMS OF M EQUATIONS   *
C     *   GIVEN THE ASYMPTOTIC EXPANSION                    *
C     *                                                     *
C     *    T(K) = SUM(A(J)*(D**(J-1))Y/FACT(J-1)*H**(J-1))  *
C     *                        J = Q+1,...,Q+P*K            *
C     *                                                     *
C     *   EVALUATED  ON  A GENERAL MESH  X(I) , I=1,...,N+1 *
C     *   D02RAX  WILL PRODUCE  S(1),...,S(N  ) = AN        *
C     *   H**(Q+P*K)  ORDER APPROXIMATION TO  T(K)    MIDWAY*
C     *   BETWEEN EACH PAIR OF CONSECUTIVE GRID POINTS.     *
C     *   FOR FIXED INTEGERS  N,P,Q, A RESTRICTION ON K IS  *
C     *                                                     *
C     *                  K .LE. (N+1-Q)/P                   *
C     *                                                     *
C     *   ALSO  P .GE. 1  AND  K .GE. 1                     *
C     *   IERROR = 1  MEANS THAT ONE OF THESE CONDITIONS HAS*
C     *   BEEN VIOLATED AND NO CORRECTION HAS BEEN COMPUTED *
C     *   A(1),...,A(Q)  ARE SET TO ZERO BY  D02RAX.        *
C     *   BOTH  Y  AND  S  ARE STORED AS VECTORS =          *
C     *   Y(1,X(1)),Y(2,X(1)),...                           *
C     *   SUBROUTINE  D02RAW  IS REQUIRED.                  *
C     *                                                     *
C     *   APRIL 1973              M. LENTINI + V. PEREYRA   *
C     *                                                     *
C     *******************************************************
C
C
C     .. Scalar Arguments ..
      INTEGER           IERROR, K, M, MTNMAX, N, NMAX, P, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  A(50), ALF(50), C(50), S(MTNMAX), X(NMAX),
     *                  Y(MTNMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACUM
      INTEGER           I, I1, IAD, II, II1, IM, IMPL, IT, ITE, ITPL, J,
     *                  JM1ETC, KK, KK1, KMID, KMID1, KMIDP1, L, MTIETC,
     *                  NF
C     .. External Subroutines ..
      EXTERNAL          D02RAW
C     .. Executable Statements ..
C
C     ***    ERROR  EXIT    ***
      IF (K.GT.(N+1-Q)/P .OR. P.LT.1 .OR. K.LT.1) GO TO 260
      IF (Q.EQ.0) GO TO 40
      DO 20 I = 1, Q
         A(I) = 0.D0
   20 CONTINUE
   40 KK1 = Q + P*K
      KK = KK1 - 1
      KMID = KK1/2
      IERROR = 0
      KMID1 = KMID - 1
C
C     UNSYMMETRIC APPROXIMATION  LEFT BOUNDARY
C
      IF (KMID1.LT.1) GO TO 120
      DO 100 I = 1, KMID1
         ITE = I
         CALL D02RAW(ITE,KK1,ITE,C,A,X,NMAX,.5D0*(X(I+1)+X(I)),ALF)
         IM = (I-1)*M
         DO 80 L = 1, M
            ACUM = 0.D0
            DO 60 J = 1, KK1
               JM1ETC = M*(J-1) + L
               ACUM = ACUM + C(J)*Y(JM1ETC)
   60       CONTINUE
            IMPL = IM + L
            S(IMPL) = ACUM
   80    CONTINUE
  100 CONTINUE
C
C     CENTER RANGE
C
  120 NF = N + 1 - KK1 + KMID
      DO 180 I = KMID, NF
         ITE = I
         CALL D02RAW(ITE,KK1,KMID,C,A,X,NMAX,.5D0*(X(I+1)+X(I)),ALF)
         I1 = I - 1
         II = I1 - KMID
         IT = I1*M
         DO 160 L = 1, M
            ACUM = 0.D0
            DO 140 J = 1, KK1
               MTIETC = M*(II+J) + L
               ACUM = ACUM + C(J)*Y(MTIETC)
  140       CONTINUE
            ITPL = IT + L
            S(ITPL) = ACUM
  160    CONTINUE
  180 CONTINUE
C
C     RIGHT  BOUNDARY
C
      KMIDP1 = KMID + 1
      II = N - KK
      II1 = II - 1
      DO 240 I = KMIDP1, KK
         IAD = II + I
         ITE = I
         CALL D02RAW(IAD,KK1,ITE,C,A,X,NMAX,.5D0*(X(IAD+1)+X(IAD)),ALF)
         IT = (II1+I)*M
         DO 220 L = 1, M
            ACUM = 0.D0
            DO 200 J = 1, KK1
               MTIETC = M*(II1+J) + L
               ACUM = ACUM + C(J)*Y(MTIETC)
  200       CONTINUE
            ITPL = IT + L
            S(ITPL) = ACUM
  220    CONTINUE
  240 CONTINUE
C
C     ...   REGULAR  EXIT  ......
C
      RETURN
  260 IERROR = 1
      RETURN
      END
