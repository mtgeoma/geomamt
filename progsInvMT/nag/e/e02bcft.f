      SUBROUTINE E02BCF(NCAP7,K,C,X,LEFT,S,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************
C     *                                                *
C     *         NAG LIBRARY SUBROUTINE  E02BCF         *
C     *                                                *
C     *  EVALUATION OF CUBIC SPLINE AND ITS            *
C     *  DERIVATIVES FROM ITS B-SPLINE REPRESENTATION  *
C     *                                                *
C     *  ROUTINE CREATED ... 17 NOV 1977               *
C     *  LATEST UPDATE ....  24 APR 1978               *
C     *  RELEASE NUMBER ...  01                        *
C     *  AUTHORS ... MAURICE G. COX AND                *
C     *              J. GEOFFREY HAYES, N.P.L.         *
C     *                                                *
C     **************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02BCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           IFAIL, LEFT, NCAP7
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NCAP7), K(NCAP7), S(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  C1, C2, C3, C4, D1N41, D1N42, D1N43, D1N44,
     *                  D2N41, D2N42, D2N43, D2N44, D3N41, D3N42, D3N43,
     *                  D3N44, E2, E3, E4, E5, HALF, K1, K2, K3, K4, K5,
     *                  K6, M11, M21, M22, M32, N41, N42, N43, N44, P4,
     *                  P5, P6, SIX, THREE
      INTEGER           IERROR, J, J1, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      HALF = 0.5D+00
      THREE = 3.0D+00
      SIX = 6.0D+00
C
C     *********  DATA VALIDATION  *********
C
C     CHECK WHETHER AT LEAST ONE INTERVAL HAS BEEN SPECIFIED -
C
      IERROR = 1
      IF (NCAP7.LT.8) GO TO 80
C
C     CHECK WHETHER THE RANGE OF DEFINITION OF THE SPLINE IS
C     STRICTLY POSITIVE IN LENGTH -
C
      IERROR = 2
      IF (K(4).GE.K(NCAP7-3)) GO TO 80
C
C     CHECK WHETHER  X  IS A VALID ARGUMENT -
C
      IF (X.LT.K(4) .OR. X.GT.K(NCAP7-3)) GO TO 80
C
C     *********  COMPUTATION  *********
C
C     BINARY SEARCH FOR INTERVAL CONTAINING  X  -
C
C     IF RIGHT-HAND DERIVATIVES ARE REQUIRED  (LEFT .NE. 1)
C     SEARCH FOR  J  SATISFYING
C     K(J + 3) .LE. X .LT. K(J + 4)
C     (SETTING  J = NCAP  IN THE EXCEPTIONAL CASE  X = K(NCAP +
C     4)).
C
C     IF LEFT-HAND DERIVATIVES ARE REQUIRED  (LEFT .EQ. 1)
C     SEARCH FOR  J  SATISFYING
C     K(J + 3) .LT. X .LE. K(J + 4)
C     (SETTING  J = 1  IN THE EXCEPTIONAL CASE  X = K(4)).
C
      IERROR = 0
      J1 = 4
      J = NCAP7 - 3
   20 L = (J1+J)/2
      IF (J-J1.LE.1) GO TO 60
      IF (LEFT.NE.1 .AND. X.GE.K(L)) GO TO 40
      IF (LEFT.EQ.1 .AND. X.GT.K(L)) GO TO 40
      J = L
      GO TO 20
   40 J1 = L
      GO TO 20
   60 J = J - 4
C
C     FORM CERTAIN CONSTANTS -
C
      K1 = K(J+1)
      K2 = K(J+2)
      K3 = K(J+3)
      K4 = K(J+4)
      K5 = K(J+5)
      K6 = K(J+6)
      E2 = X - K2
      E3 = X - K3
      E4 = K4 - X
      E5 = K5 - X
      P4 = K4 - K1
      P5 = K5 - K2
      P6 = K6 - K3
C
C     FORM BASIS FUNCTIONS AND THEIR DERIVATIVES -
C
C     THE VALUES OF THE NON-ZERO UN-NORMALIZED B-SPLINES OF ORDER R
C     ARE DENOTED BY  MR1, MR2,...  .  THE CORRESPONDING NORMALIZED
C     B-SPLINES ARE DENOTED BY  NR1, NR2,...  AND THEIR DERIVATIVES
C     OF ORDER  L  BY  DLNR1, DLNR2,...  .
C
      M11 = SIX/(K4-K3)
      M21 = -M11/(K4-K2)
      M22 = M11/(K5-K3)
      D3N41 = M21/P4
      M32 = (M21-M22)/P5
      D3N44 = M22/P6
      D3N42 = -D3N41 - M32
      D3N43 = M32 - D3N44
      M21 = -E4*M21
      M22 = E3*M22
      D2N41 = M21/P4
      M32 = (M21-M22)/P5
      D2N44 = M22/P6
      D2N42 = -D2N41 - M32
      D2N43 = M32 - D2N44
      M21 = HALF*M21
      M22 = HALF*M22
      D1N41 = -E4*M21/P4
      M32 = (E2*M21+E5*M22)/P5
      D1N44 = E3*M22/P6
      D1N42 = -D1N41 - M32
      D1N43 = M32 - D1N44
      N41 = -E4*D1N41/THREE
      N42 = (-(X-K1)*D1N41+E5*M32)/THREE
      N43 = (E2*M32+(K6-X)*D1N44)/THREE
      N44 = E3*D1N44/THREE
C
C     FORM THE VALUES OF THE CUBIC SPLINE AND ITS DERIVATIVES -
C
      C1 = C(J)
      C2 = C(J+1)
      C3 = C(J+2)
      C4 = C(J+3)
      S(1) = C1*N41 + C2*N42 + C3*N43 + C4*N44
      S(2) = C1*D1N41 + C2*D1N42 + C3*D1N43 + C4*D1N44
      S(3) = C1*D2N41 + C2*D2N42 + C3*D2N43 + C4*D2N44
      S(4) = C1*D3N41 + C2*D3N42 + C3*D3N43 + C4*D3N44
C
C     *********  ERROR DIAGNOSTICS  *********
C
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
C
C     END OF SUBROUTINE  E02BCF
C
      END
