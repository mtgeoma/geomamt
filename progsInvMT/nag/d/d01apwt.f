      SUBROUTINE D01APW(ALFA,BETA,RI,RJ,RG,RH,INTEGR)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QMOMO.
C     ..................................................................
C
C        PURPOSE
C           THIS ROUTINE COMPUTES MODIFIED CHEBYSHEV MOMENTS.
C           THE K-TH MODIFIED CHEBYSHEV MOMENT IS DEFINED AS THE
C           INTEGRAL OVER (-1,1) OF W(X)*T(K,X), WHERE T(K,X) IS THE
C           CHEBYSHEV POLYNOMIAL OF DEGREE K.
C
C        PARAMETERS
C           ALFA   - REAL
C                    PARAMETER IN THE WEIGHT FUNCTION W(X), ALFA.GT.(-1)
C
C           BETA   - REAL
C                    PARAMETER IN THE WEIGHT FUNCTION W(X), BETA.GT.(-1)
C
C           RI     - REAL
C                    VECTOR OF DIMENSION 25
C                    RI(K) IS THE INTEGRAL OVER (-1,1) OF
C                    (1+X)**ALFA*T(K-1,X), K = 1, ..., 25.
C
C           RJ     - REAL
C                    VECTOR OF DIMENSION 25
C                    RJ(K) IS THE INTEGRAL OVER (-1,1) OF
C                    (1-X)**BETA*T(K-1,X), K = 1, ..., 25.
C
C           RG     - REAL
C                    VECTOR OF DIMENSION 25
C                    RG(K) IS THE INTEGRAL OVER (-1,1) OF
C                    (1+X)**ALFA*LOG((1+X)/2)*T(K-1,X), K = 1, ...,25.
C
C           RH     - REAL
C                    VECTOR OF DIMENSION 25
C                    RH(K) IS THE INTEGRAL OVER (-1,1) OF
C                    (1-X)**BETA*LOG((1-X)/2)*T(K-1,X), K = 1, ..., 25.
C
C           INTEGR - INTEGER
C                    INPUT PARAMETER INDICATING THE MODIFIED MOMENTS
C                    TO BE COMPUTED
C                    INTEGR = 1 COMPUTE RI, RJ
C                           = 2 COMPUTE RI, RJ, RG
C                           = 3 COMPUTE RI, RJ, RH
C                           = 4 COMPUTE RI, RJ, RG, RH
C
C     ..................................................................
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, BETA
      INTEGER           INTEGR
C     .. Array Arguments ..
      DOUBLE PRECISION  RG(25), RH(25), RI(25), RJ(25)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFP1, ALFP2, AN, ANM1, BETP1, BETP2, RALF, RBET
      INTEGER           I, IM1
C     .. Executable Statements ..
      ALFP1 = ALFA + 1.0D+00
      BETP1 = BETA + 1.0D+00
      ALFP2 = ALFA + 2.0D+00
      BETP2 = BETA + 2.0D+00
      RALF = 2.0D+00**ALFP1
      RBET = 2.0D+00**BETP1
C
C           COMPUTE RI, RJ USING A FORWARD RECURRENCE RELATION.
C
      RI(1) = RALF/ALFP1
      RJ(1) = RBET/BETP1
      RI(2) = RI(1)*ALFA/ALFP2
      RJ(2) = RJ(1)*BETA/BETP2
      AN = 2.0D+00
      ANM1 = 1.0D+00
      DO 20 I = 3, 25
         RI(I) = -(RALF+AN*(AN-ALFP2)*RI(I-1))/(ANM1*(AN+ALFP1))
         RJ(I) = -(RBET+AN*(AN-BETP2)*RJ(I-1))/(ANM1*(AN+BETP1))
         ANM1 = AN
         AN = AN + 1.0D+00
   20 CONTINUE
      IF (INTEGR.EQ.1) GO TO 120
      IF (INTEGR.EQ.3) GO TO 60
C
C           COMPUTE RG USING A FORWARD RECURRENCE RELATION.
C
      RG(1) = -RI(1)/ALFP1
      RG(2) = -(RALF+RALF)/(ALFP2*ALFP2) - RG(1)
      AN = 2.0D+00
      ANM1 = 1.0D+00
      IM1 = 2
      DO 40 I = 3, 25
         RG(I) = -(AN*(AN-ALFP2)*RG(IM1)-AN*RI(IM1)+ANM1*RI(I))
     *           /(ANM1*(AN+ALFP1))
         ANM1 = AN
         AN = AN + 1.0D+00
         IM1 = I
   40 CONTINUE
      IF (INTEGR.EQ.2) GO TO 120
C
C           COMPUTE RH USING A FORWARD RECURRENCE RELATION.
C
   60 RH(1) = -RJ(1)/BETP1
      RH(2) = -(RBET+RBET)/(BETP2*BETP2) - RH(1)
      AN = 2.0D+00
      ANM1 = 1.0D+00
      IM1 = 2
      DO 80 I = 3, 25
         RH(I) = -(AN*(AN-BETP2)*RH(IM1)-AN*RJ(IM1)+ANM1*RJ(I))
     *           /(ANM1*(AN+BETP1))
         ANM1 = AN
         AN = AN + 1.0D+00
         IM1 = I
   80 CONTINUE
      DO 100 I = 2, 25, 2
         RH(I) = -RH(I)
  100 CONTINUE
  120 DO 140 I = 2, 25, 2
         RJ(I) = -RJ(I)
  140 CONTINUE
      RETURN
      END
