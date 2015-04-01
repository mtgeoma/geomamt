      SUBROUTINE D02NNW(NEQ,T,TOUT,H0,Y,YDOTI,EWT,RTOL,ATOL,ITOL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     OLD NAME ITSTEP
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H0, T, TOUT
      INTEGER           ITOL, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), EWT(NEQ), RTOL(*), Y(NEQ), YDOTI(NEQ)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           IDEV, INORM, IOVFLO, ITRACE
C     .. Local Scalars ..
      DOUBLE PRECISION  ATOLI, AYI, SUM, TDIST, TOL, W0
      INTEGER           I, IFZAF
C     .. External Functions ..
      DOUBLE PRECISION  D02ZAF
      EXTERNAL          D02ZAF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, SQRT
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /FD02NM/DUNFLO, UROUND, IOVFLO
      COMMON            /HD02NM/INORM
C     .. Save statement ..
      SAVE              /HD02NM/, /FD02NM/, /AD02NM/
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C THE CODING BELOW COMPUTES THE STEP SIZE, H0, TO BE ATTEMPTED ON THE
C FIRST STEP, UNLESS THE USER HAS SUPPLIED A VALUE FOR THIS.
C FIRST CHECK THAT TOUT - T DIFFERS SIGNIFICANTLY FROM ZERO.
C A SCALAR TOLERANCE QUANTITY TOL IS COMPUTED, AS MAX(RTOL(I))
C IF THIS IS POSITIVE, OR MAX(ATOL(I)/ABS(Y(I))) OTHERWISE, ADJUSTED
C SO AS TO BE BETWEEN 100*UROUND AND 1.0E-3.
C THEN THE COMPUTED VALUE H0 IS GIVEN BY..
C
C   H0**2 = TOL / ( W0**-2 + '' YDOT / YWT '' **2  )
C
C WHERE   W0      = MAX ( ABS(T), ABS(TOUT) ),
C         YDOT(I) = I-TH COMPONENT OF INITIAL VALUE OF DY/DT,
C         YWT(I)  = EWT(I)/TOL  (A WEIGHT FOR Y(I)).
C  AND    THE NORM USED  '' . '' IS THAT CHOSEN BY THE USER.
C THE SIGN OF H0 IS INFERRED FROM THE INITIAL VALUES OF TOUT AND T.
C-----------------------------------------------------------------------
      IF (H0.NE.0.0D0) GO TO 100
      TDIST = ABS(TOUT-T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST.LT.2.0D0*UROUND*W0) THEN
         CALL D02NNQ(
     *' TOUT (=R1) TOO CLOSE TO T (=R2) TO                  START INTEGR
     *ATION',1,0,0,0,2,TOUT,T)
         GO TO 100
      END IF
      TOL = RTOL(1)
      IF (ITOL.LE.2) GO TO 40
      DO 20 I = 1, NEQ
         TOL = MAX(TOL,RTOL(I))
   20 CONTINUE
   40 IF (TOL.GT.0.0D0) GO TO 80
      ATOLI = ATOL(1)
      DO 60 I = 1, NEQ
         IF (ITOL.EQ.2 .OR. ITOL.EQ.4) ATOLI = ATOL(I)
         AYI = ABS(Y(I))
         IF (AYI.NE.0.0D0) TOL = MAX(TOL,ATOLI/AYI)
   60 CONTINUE
   80 TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      IFZAF = 1
      SUM = D02ZAF(NEQ,YDOTI,EWT,IFZAF)
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
  100 CONTINUE
      RETURN
      END
