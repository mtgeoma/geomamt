      SUBROUTINE F04AMF(A,IA,X,IX,B,IB,M,N,IP,ETA,QR,IQR,ALPHA,E,Y,Z,R,
     *                  IPIVOT,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16A REVISED. IER-1003 (JUN 1993).
C     LEAST SQUARES SOLUTION
C     THE ARRAY A(M,N) CONTAINS THE GIVEN MATRIX OF AN
C     OVERDETERMINED SYSTEM OF M LINEAR EQUATIONS IN N UNKNOWNS
C     (M.GE.N). FOR THE IP RIGHT HAND SIDES GIVEN AS THE COLUMNS
C     OF THE ARRAY B(M,IP), THE LEAST SQUARES SOLUTIONS ARE
C     COMPUTED AND STORED AS THE COLUMNS OF THE ARRAY X(N,IP).
C     IF RANK(A).LT.N THEN THE PROBLEM IS LEFT UNSOLVED AND IFAIL
C     IS SET EQUAL TO 1. IN EITHER CASE A AND B ARE
C     LEFT INTACT. ETA IS THE RELATIVE MACHINE PRECISION.
C     ADDITIONAL PRECISION INNERPRODUCTS ARE ABSOLUTELY NECESSARY.
C     1ST. MARCH  1972
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AMF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ETA
      INTEGER           IA, IB, IFAIL, IP, IQR, IX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), ALPHA(N), B(IB,IP), E(N), QR(IQR,N),
     *                  R(M), X(IX,IP), Y(N), Z(N)
      INTEGER           IPIVOT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, D2, ETA2, NORME0, NORME1, NORMY0
      INTEGER           I, IFAIL1, ISAVE, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01AXF, F04ANF, X03AAF
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      IFAIL1 = 0
      DO 40 J = 1, N
         DO 20 I = 1, M
            QR(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
      CALL F01AXF(M,N,QR,IQR,ALPHA,IPIVOT,Y,E,IFAIL)
      IF (IFAIL.EQ.0) GO TO 60
      IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
   60 ETA2 = ETA*ETA
      DO 260 K = 1, IP
C        SOLUTION FOR K-TH RIGHT HAND SIDE
         DO 80 I = 1, M
            R(I) = B(I,K)
   80    CONTINUE
         CALL F04ANF(M,N,QR,IQR,ALPHA,IPIVOT,R,Y,Z)
         DO 100 I = 1, M
            CALL X03AAF(A(I,1),N*IA-I+1,Y(1),N,N,IA,1,-B(I,K)
     *                  ,0.0D0,D1,D2,.TRUE.,IFAIL1)
            R(I) = -D1
  100    CONTINUE
         CALL F04ANF(M,N,QR,IQR,ALPHA,IPIVOT,R,E,Z)
         NORMY0 = DNRM2(N,Y,1)
         NORME1 = DNRM2(N,E,1)
         IF (NORME1.LE.0.25D0*NORMY0) GO TO 140
         IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
         RETURN
C        NO ATTEMPT AT OBTAINING THE SOLUTION IS MADE UNLESS THE NORM
C        OF THE FIRST CORRECTION IS SIGNIFICANTLY SMALLER THAN THE
C        NORM OF THE INITIAL SOLUTION
C        LABEL IMPROVE
  140    DO 160 I = 1, N
            Y(I) = Y(I) + E(I)
  160    CONTINUE
         IF (NORME1.LE.ETA*NORMY0) GO TO 220
C        TERMINATE THE ITERATION IF THE
C        CORRECTION WAS OF LITTLE SIGNIFICANCE
         DO 180 I = 1, M
            CALL X03AAF(A(I,1),N*IA-I+1,Y(1),N,N,IA,1,-B(I,K)
     *                  ,0.0D0,D1,D2,.TRUE.,IFAIL1)
            R(I) = -D1
  180    CONTINUE
         CALL F04ANF(M,N,QR,IQR,ALPHA,IPIVOT,R,E,Z)
         NORME0 = NORME1
         NORME1 = DNRM2(N,E,1)
         IF (NORME1.LE.0.25D0*NORME0) GO TO 140
C        TERMINATE THE ITERATION ALSO IF THE NORM OF THE CORRECTION
C        FAILED TO DECREASE SUFFICIENTLY AS COMPARED WITH THE NORM
C        OF THE PREVIOUS CORRECTION
C        LABEL STORE
  220    DO 240 I = 1, N
            X(I,K) = Y(I)
  240    CONTINUE
C        END OF THE K-TH RIGHT HAND SIDE
  260 CONTINUE
      RETURN
      END
