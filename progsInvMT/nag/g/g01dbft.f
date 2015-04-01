      SUBROUTINE G01DBF(N,PP,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     AS APPL. STATIST. ALGORITHM AS 177.3 (1982), VOL. 31
C
C     COMPUTES AN APPROXIMATION TO THE EXPECTED VALUES OF NORMAL ORDER
C     STATISTICS FOR A SAMPLE OF SIZE N.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01DBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  PP(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, AL1, AN, B1, BB, D, E1, E2, ZERO
      INTEGER           I, IERROR, J, K, M, N2
C     .. Local Arrays ..
      DOUBLE PRECISION  ALAM(4), DL1(4), DL2(4), EPS(4), GAM(4)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01CEF, G01DBZ
      INTEGER           P01ABF
      EXTERNAL          G01CEF, G01DBZ, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Data statements ..
C
C      INITIALISE CONSTANTS
C
      DATA              EPS(1), EPS(2), EPS(3), EPS(4)/0.419885D0,
     *                  0.450536D0, 0.456936D0, 0.468488D0/, DL1(1),
     *                  DL1(2), DL1(3), DL1(4)/0.112063D0, 0.121770D0,
     *                  0.239299D0, 0.215159D0/, DL2(1), DL2(2), DL2(3),
     *                  DL2(4)/0.080122D0, 0.111348D0, -0.211867D0,
     *                  -0.115049D0/, GAM(1), GAM(2), GAM(3),
     *                  GAM(4)/0.474798D0, 0.469051D0, 0.208597D0,
     *                  0.259784D0/, ALAM(1), ALAM(2), ALAM(3),
     *                  ALAM(4)/0.282765D0, 0.304856D0, 0.407708D0,
     *                  0.414093D0/, BB/-0.283833D0/, D/-0.106136D0/,
     *                  B1/0.5641896D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      IERROR = 0
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) '** ON ENTRY, THE VALUE OF N (', N,
     *     ' ) IS LESS THAN ONE'
         GO TO 100
      END IF
      IF (N.EQ.1) THEN
         PP(1) = 0.0D0
         GO TO 80
      END IF
      N2 = N/2
      M = N2
      IF (2*N2.NE.N) M = M + 1
      PP(1) = -B1
      PP(2) = B1
      IF (N.EQ.2) GO TO 80
C
C     CALCULATE NORMAL AREAS FOR 3 LARGEST RANKITS
C
      AN = N
      K = MIN(3,N2)
      DO 20 I = 1, K
         AI = I
         E1 = (AI-EPS(I))/(AN+GAM(I))
         E2 = E1**ALAM(I)
         J = N + 1 - I
         PP(J) = E1 + E2*(DL1(I)+E2*DL2(I))/AN - G01DBZ(I,N)
   20 CONTINUE
C
C     CALCULATE AREAS FOR REMAINING RANKITS
C
      DO 40 I = 4, N2
         AI = I
         AL1 = ALAM(4) + BB/(AI+D)
         E1 = (AI-EPS(4))/(AN+GAM(4))
         E2 = E1**AL1
         J = N + 1 - I
         PP(J) = E1 + E2*(DL1(4)+E2*DL2(4))/AN - G01DBZ(I,N)
   40 CONTINUE
C
C     CONVERT NORMAL TAIL AREAS TO NORMAL DEVIATES
C
      IF (M.NE.N2) PP(M) = ZERO
      DO 60 I = 1, N2
         J = M + I
         PP(J) = -G01CEF(PP(J),IERROR)
         K = N2 - I + 1
         PP(K) = -PP(J)
   60 CONTINUE
   80 IFAIL = 0
      RETURN
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
      RETURN
C
99999 FORMAT (1X,A,I16,A)
      END
