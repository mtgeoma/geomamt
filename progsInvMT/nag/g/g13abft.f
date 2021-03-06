      SUBROUTINE G13ABF(X,NX,NK,XM,XV,R,STAT,IFAIL)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-512 (AUG 1986).
C
C     G13ABF COMPUTES THE SAMPLE AUTOCORRELATION FUNCTION
C     OF AN INPUT TIME SERIES. IT ALSO COMPUTES THE SAMPLE MEAN,
C     THE SAMPLE VARIANCE, AND A STATISTIC WHICH MAY BE USED TO
C     TEST THE HYPOTHESIS THAT THE TRUE AUTOCORRELATION
C     FUNCTION IS ZERO.
C
C     CONTRIBUTORS - G. TUNNICLIFFE WILSON, C. DALY (LANCASTER
C     UNIV.)
C     VALIDATOR    - T. LAMBERT (NAG CENTRAL OFFICE)
C
C     USES NAG LIBRARY ROUTINE P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13ABF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STAT, XM, XV
      INTEGER           IFAIL, NK, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NK), X(NX)
C     .. Local Scalars ..
      DOUBLE PRECISION  RL, RNX, SSX, XDIFF, ZERO, T2, T3
      INTEGER           I, IERROR, IH, J, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          P01ABF, X02AJF
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
C
C     TEST FOR ERRORS IN INPUT PARAMETERS
C
      IERROR = 1
      IF (NK.LE.0) GO TO 100
      IF (NX.LE.NK) GO TO 100
C
C     CALCULATE MEAN AND SAMPLE VARIANCE
C
      RNX = NX
      XV = ZERO
      STAT = ZERO
      XM = X(1)
      SSX = ZERO
      DO 40 I = 2, NX
C
         T2 = X(I) - XM
         T3 = T2/I
         XM = XM + T3
         SSX = SSX + (I-1)*T2*T3
C
   40 CONTINUE
C
C     TEST FOR ZERO VARIANCE
C
      IERROR = 2
      XV = SSX/(RNX-1.0D0)
      IF (XV.LE.4.0D0*X02AJF()**2) GO TO 100
C
C     CALCULATE SAMPLE AUTOCORRELATION COEFFICIENTS
C     AND TEST STATISTIC STAT
C
      DO 80 L = 1, NK
         RL = ZERO
         IH = NX - L
         DO 60 I = 1, IH
            J = I + L
            RL = RL + (X(I)-XM)*(X(J)-XM)
   60    CONTINUE
         RL = RL/SSX
         STAT = STAT + RL*RL*RNX
         R(L) = RL
   80 CONTINUE
      IFAIL = 0
      RETURN
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
