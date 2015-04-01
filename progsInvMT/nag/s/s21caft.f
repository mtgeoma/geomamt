      SUBROUTINE S21CAF(U,M,SN,CN,DN,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Returns Jacobian elliptic functions sn, cn, dn.
C     Based on routine by R. Bulirsch (Num. Math. 7, pp 76-90, 1965).
C
C     .. Parameters ..
      INTEGER           MAXIT
      PARAMETER         (MAXIT=30)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='S21CAF')
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CN, DN, M, SN, U
      INTEGER           IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ABSU, B, C, D, EPS, LARGE, MC, ROOTLG, SMALL,
     *                  X
      INTEGER           I, IERR, K, NREC
      LOGICAL           NEGATE
C     .. Local Arrays ..
      DOUBLE PRECISION  MU(0:MAXIT), NU(0:MAXIT)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, COS, COSH, LOG, SIN, SQRT, TANH
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      SMALL = X02AMF()
      LARGE = ONE/SMALL
      ROOTLG = SQRT(LARGE)
      ABSU = ABS(U)
      IF (ABSU.GT.ROOTLG) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99999) ABSU, ROOTLG
         SN = ZERO
         CN = ZERO
         DN = ZERO
      ELSE IF (ABSU.LT.ONE/ROOTLG) THEN
         IF (ABS(M).GT.ROOTLG) THEN
            IERR = 2
            NREC = 2
            WRITE (REC,FMT=99998) ABS(M), ROOTLG
            SN = ZERO
            CN = ZERO
            DN = ZERO
         ELSE
C           Use leading terms of power series (Abramowitz and Stegun,
C           16.22, page 575).
            SN = U
            CN = ONE
            DN = ONE
         END IF
      ELSE IF (M.EQ.ZERO) THEN
         SN = SIN(U)
         CN = COS(U)
         DN = ONE
      ELSE IF (M.EQ.ONE) THEN
         SN = TANH(U)
         IF (ABSU.LT.LOG(LARGE)) THEN
            CN = ONE/COSH(U)
         ELSE
C           Avoids overflow in call to COSH.
            CN = ZERO
         END IF
         DN = CN
      ELSE
         MC = ONE - M
         EPS = SQRT(X02AJF())
         X = U
C
         NEGATE = MC .LT. ZERO
         IF (NEGATE) THEN
C           Change parameter (Abramowitz and Stegun, 16.10, page 573)
            MC = -MC/M
            D = SQRT(M)
            X = D*X
         END IF
C
         A = ONE
         I = 0
   20    MU(I) = A
         MC = SQRT(MC)
         NU(I) = MC
         C = HALF*(A+MC)
         IF (ABS(A-MC).GT.EPS*A) THEN
            MC = A*MC
            A = C
            IF (I.LT.MAXIT) THEN
               I = I + 1
               GO TO 20
            END IF
         END IF
C
         X = C*X
         SN = SIN(X)
         CN = COS(X)
         DN = ONE
C
         IF (SN.NE.ZERO) THEN
            A = CN/SN
            C = A*C
            DO 40 K = I, 0, -1
               B = MU(K)
               A = C*A
               C = DN*C
               DN = (NU(K)+A)/(B+A)
               A = C/B
   40       CONTINUE
            A = ONE/SQRT(C*C+ONE)
            IF (SN.LT.ZERO) THEN
               SN = -A
            ELSE
               SN = A
            END IF
            CN = C*SN
         END IF
C
         IF (NEGATE) THEN
            A = DN
            DN = CN
            CN = A
            SN = SN/D
         END IF
C
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, abs(U) is too large : abs(U) = ',D13.5,
     *       /'    it must be less than ',D13.5)
99998 FORMAT (' ** On entry, abs(M) is too large when used in conjunct',
     *       'ion with the supplied',/'    argument U : abs(M) = ',
     *       D13.5,' it must be less than',D13.5)
      END
