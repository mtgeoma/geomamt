      SUBROUTINE E02AGF(M,KPLUS1,NROWS,XMIN,XMAX,X,Y,W,MF,XF,YF,LYF,IP,
     *                  A,S,NP1,WRK,LWRK,IWRK,LIWRK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-314 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14C REVISED. IER-878 (NOV 1990).
C
C     **************************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE CONFIT
C
C     CREATED  20/8/1979      UPDATED  16/5/80     RELEASE NO. 00/05
C
C     AUTHOR... GERALD T ANTHONY.
C     NATIONAL PHYSICAL LABORATORY,
C     TEDDINGTON, MIDDLESEX TW11 0LW,
C     ENGLAND.
C
C     **************************************************************
C
C     E02AGF CALLS E01AEW TO DETERMINE A POLYNOMIAL MU(X) WHICH
C     INTERPOLATES THE GIVEN CONSTRAINTS AND A POLYNOMIAL NU(X)
C     WHICH HAS VALUE ZERO WHERE A CONSTRAINED VALUE IS SPECIFIED.
C     IT THEN CALLS E02ADZ TO FIT Y-MU(X) AS A POLYNOMIAL IN X
C     WITH FACTOR NU(X). FINALLY THE COEFFICIENTS OF MU ARE ADDED
C     TO THESE FITS TO GIVE THE COEFFICIENTS OF THE CONSTRAINED
C     FITS TO Y. ALL POLYNOMIALS ARE EXPRESSED IN CHEBYSHEV
C     SERIES FORM
C
C     INPUT PARAMETERS
C        M        THE NUMBER OF DATA POINTS TO BE FITTED
C        KPLUS1   FITS WITH UP TO KPLUS1 COEFFICIENTS ARE REQUIRED
C        NROWS    FIRST DIMENSION OF ARRAY A WHERE COEFFICIENTS ARE
C                 TO BE STORED
C        XMIN,    END POINTS OF THE RANGE OF THE
C        XMAX     INDEPENDENT VARIABLE
C        X, Y, W  ARRAYS OF DATA VALUES OF THE INDEPENDENT VARIABLE,
C                 DEPENDENT VARIABLE AND WEIGHT, RESPECTIVELY
C        MF       NUMBER OF X VALUES AT WHICH A CONSTRAINT
C                 IS SPECIFIED
C        XF       ARRAY OF VALUES OF THE INDEPENDENT
C                 VARIABLE AT WHICH CONSTRAINTS ARE
C                 SPECIFIED
C        YF       ARRAY OF SPECIFIED VALUES AND DERIVATIVES OF THE
C                 DEPENDENT VARIABLE IN THE ORDER
C                    Y1, Y1 DERIV, Y1 2ND DERIV,...., Y2,....
C        LYF      DIMENSION OF ARRAY YF
C        IP       INTEGER ARRAY OF DEGREES OF DERIVATIVES
C                 SPECIFIED AT EACH POINT XF
C
C     OUTPUT PARAMETERS
C        A        ON EXIT, 2 PARAMETER ARRAY CONTAINING THE
C                 COEFFICIENTS OF THE CHEBYSHEV SERIES
C                 REPRESENTATION OF THE FITS, A(I+1, J+1)
C                 CONTAINS THE COEFFICIENT OF TJ IN THE FIT
C                 OF DEGREE I, I = N,N+1,...,K,  J =
C                 0,1,...,I  WHERE N = NP1 - 1
C        S        ON EXIT, ARRAY CONTAINING THE R.M.S. RESIDUAL FOR
C                 EACH DEGREE OF FIT FROM N TO K
C        NP1      ON EXIT, CONTAINS N + 1, WHERE N IS THE
C                 TOTAL NUMBER OF INTERPOLATION CONDITIONS
C
C        IFAIL    FAILURE INDICATOR
C                  0 - SUCCESSFUL TERMINATION
C                  1 - AT LEAST ONE OF THE FOLLOWING CONDITIONS
C                      HAS BEEN VIOLATED
C                      LYF    AT LEAST N
C                      LWRK   AT LEAST 2*N + 2 + THE LARGER OF
C                             4*M + 3*KPLUS1 AND 8*NP1 +
C                             5*IMAX + MF - 3 WHERE IMAX =
C                             1 + MAX(IP(I))
C                      LIWRK  AT LEAST 2*MF + 2
C                      KPLUS1 AT LEAST NP1
C                      M      AT LEAST 1
C                      NROWS  AT LEAST KPLUS1
C                      MF     AT LEAST 1
C                  2 - FOR SOME I, IP(I) IS LESS THAN 0
C                  3 - XMIN IS NOT STRICTLY LESS THAN XMAX
C                      OR FOR SOME I, XF(I) IS NOT IN RANGE
C                      XMIN TO XMAX OR THE XF(I) ARE NOT
C                      DISTINCT
C                  4 - FOR SOME I, X(I) IS NOT IN RANGE XMIN TO XMAX
C                  5 - THE X(I) ARE NOT NON-DECREASING
C                  6 - THE NUMBER OF DISTINCT VALUES OF X(I) WITH
C                      NON-ZERO WEIGHT IS LESS THAN KPLUS1 - NP1
C                  7 - E01AEW HAS FAILED TO CONVERGE, IE
C                      THE CONSTRAINT CANNOT BE SATISFIED
C                      WITH SUFFICIENT ACCURACY
C
C     WORKSPACE PARAMETERS
C        WRK      REAL WORKSPACE ARRAY
C        LWRK     DIMENSION OF WRK.   LWRK MUST BE AT LEAST
C                 2*N + 2 + THE LARGER OF
C                 4*M + 3*KPLUS1 AND 8*NP1 + 5*IMAX + MF - 3
C                 WHERE IMAX = 1 + MAX(IP(I))
C        IWRK     INTEGER WORKSPACE ARRAY
C        LIWRK    DIMENSION OF IWRK.   LIWRK MUST BE AT LEAST
C                 2*MF + 2
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02AGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           IFAIL, KPLUS1, LIWRK, LWRK, LYF, M, MF, NP1,
     *                  NROWS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWS,KPLUS1), S(KPLUS1), W(M), WRK(LWRK),
     *                  X(M), XF(MF), Y(M), YF(LYF)
      INTEGER           IP(MF), IWRK(LIWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMUJ, ONE, XI, XMU, ZERO
      INTEGER           I, IERROR, IM1, IMAX, IPI, IYMUX, J, LW, MDIST,
     *                  N, NANU, NEPS, NSER, NWRK, NWRK1, NWRK2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01AEW, E02ADZ, E02AKF
C     .. Data statements ..
      DATA              ONE, ZERO/1.0D+0, 0.0D+0/
C     .. Executable Statements ..
      IERROR = 1
      IF (MF.LT.1) GO TO 280
      IERROR = 2
      IMAX = 0
      NP1 = 1
      DO 20 I = 1, MF
         IPI = IP(I) + 1
         IF (IPI.LT.1) GO TO 280
         IF (IPI.GT.IMAX) IMAX = IPI
         NP1 = NP1 + IPI
   20 CONTINUE
      N = NP1 - 1
      IERROR = 1
      IF (LYF.LT.N .OR. LIWRK.LT.2*MF+2) GO TO 280
      I = 4*M + 3*KPLUS1
      LW = 8*NP1 + 5*IMAX + MF - 3
      IF (LW.LT.I) LW = I
      LW = LW + 2*NP1
      IF (LW.GT.LWRK .OR. NP1.GT.KPLUS1 .OR. M.LT.1) GO TO 280
      NANU = NP1 + 1
      NWRK = NANU + NP1
      IERROR = 3
      IF (XMAX.LE.XMIN) GO TO 280
      IF (XF(1).GT.XMAX .OR. XF(1).LT.XMIN) GO TO 280
      IF (MF.EQ.1) GO TO 80
      DO 60 I = 2, MF
         XI = XF(I)
         IF (XI.GT.XMAX .OR. XI.LT.XMIN) GO TO 280
         IM1 = I - 1
         DO 40 J = 1, IM1
            IF (XI.EQ.XF(J)) GO TO 280
   40    CONTINUE
   60 CONTINUE
   80 IERROR = 6
      XMU = XMIN
      IF (X(1).EQ.XMIN) XMU = XMAX
      MDIST = 0
      DO 140 I = 1, M
         XI = X(I)
         IF (XI.EQ.XMU .OR. W(I).EQ.ZERO) GO TO 140
         DO 100 J = 1, MF
            IF (XI.EQ.XF(J)) GO TO 120
  100    CONTINUE
         MDIST = MDIST + 1
  120    XMU = XI
  140 CONTINUE
      IF (MDIST.LT.KPLUS1-N) GO TO 280
      IWRK(1) = 1
      IERROR = 1
      CALL E01AEW(MF,XMIN,XMAX,XF,YF,IP,N,NP1,5,20,WRK,WRK(NANU),
     *            WRK(NWRK),LWRK-NWRK+1,IWRK,LIWRK,IERROR)
      IF (IERROR.EQ.0) GO TO 160
      IERROR = 7
      GO TO 280
  160 DO 200 I = 1, M
         IERROR = 1
         CALL E02AKF(N,XMIN,XMAX,WRK,1,NP1,X(I),XMU,IERROR)
         IF (IERROR.EQ.0) GO TO 180
         IERROR = 4
         GO TO 280
  180    IYMUX = NWRK + I - 1
C
C        STORE Y - MU(X) AT ITH DATA POINT
C
         WRK(IYMUX) = Y(I) - XMU
  200 CONTINUE
      NWRK1 = NWRK + M
      NWRK2 = NWRK1 + 2*M
      NSER = NWRK2 + 2*KPLUS1
      NEPS = NSER + KPLUS1
      IERROR = 1
      CALL E02ADZ(1,M,M,KPLUS1,NROWS,1,1,X,WRK(NWRK),W,XMIN,XMAX,NP1,
     *            WRK(NANU),WRK(NWRK1),WRK(NWRK2),A,S,WRK(NSER),
     *            WRK(NEPS),IERROR)
      IF (IERROR.EQ.0) GO TO 220
      IERROR = IERROR + 3
      IF (IERROR.EQ.8) IERROR = 1
      GO TO 280
  220 DO 260 J = 1, N
         AMUJ = WRK(J)
         DO 240 I = NP1, KPLUS1
            A(I,J) = A(I,J) + AMUJ
  240    CONTINUE
  260 CONTINUE
  280 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
C     END E02AGF
      END
