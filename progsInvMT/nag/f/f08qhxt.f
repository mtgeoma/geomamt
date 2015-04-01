      SUBROUTINE F08QHX(LTRANS,NA,NW,SMIN,CA,A,LDA,D1,D2,B,LDB,WR,WI,X,
     *                  LDX,SCALE,XNORM,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLALN2(LTRANS,NA,NW,SMIN,CA,A,LDA,D1,D2,B,LDB,
C    *                  WR,WI,X,LDX,SCALE,XNORM,INFO)
C
C  Purpose
C  =======
C
C  DLALN2 solves a system of the form  (ca A - w D ) X = s B
C  or (ca A' - w D) X = s B   with possible scaling ("s") and
C  perturbation of A.  (A' means A-transpose.)
C
C  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
C  real diagonal matrix, w is a real or complex value, and X and B are
C  NA x 1 matrices -- real if w is real, complex if w is complex.  NA
C  may be 1 or 2.
C
C  If w is complex, X and B are represented as NA x 2 matrices,
C  the first column of each being the real part and the second
C  being the imaginary part.
C
C  "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
C  so chosen that X can be computed without overflow.  X is further
C  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
C  than overflow.
C
C  If both singular values of (ca A - w D) are less than SMIN,
C  SMIN*identity will be used instead of (ca A - w D).  If only one
C  singular value is less than SMIN, one element of (ca A - w D) will be
C  perturbed enough to make the smallest singular value roughly SMIN.
C  If both singular values are at least SMIN, (ca A - w D) will not be
C  perturbed.  In any case, the perturbation will be at most some small
C  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
C  are computed by infinity-norm approximations, and thus will only be
C  correct to a factor of 2 or so.
C
C  Note: all input quantities are assumed to be smaller than overflow
C  by a reasonable factor.  (See BIGNUM.)
C
C  Arguments
C  ==========
C
C  LTRANS  (input) LOGICAL
C          =.TRUE.:  A-transpose will be used.
C          =.FALSE.: A will be used (not transposed.)
C
C  NA      (input) INTEGER
C          The size of the matrix A.  It may (only) be 1 or 2.
C
C  NW      (input) INTEGER
C          1 if "w" is real, 2 if "w" is complex.  It may only be 1
C          or 2.
C
C  SMIN    (input) DOUBLE PRECISION
C          The desired lower bound on the singular values of A.  This
C          should be a safe distance away from underflow or overflow,
C          say, between (underflow/machine precision) and  (machine
C          precision * overflow ).  (See BIGNUM and ULP.)
C
C  CA      (input) DOUBLE PRECISION
C          The coefficient c, which A is multiplied by.
C
C  A       (input) DOUBLE PRECISION array, dimension (LDA,NA)
C          The NA x NA matrix A.
C
C  LDA     (input) INTEGER
C          The leading dimension of A.  It must be at least NA.
C
C  D1      (input) DOUBLE PRECISION
C          The 1,1 entry in the diagonal matrix D.
C
C  D2      (input) DOUBLE PRECISION
C          The 2,2 entry in the diagonal matrix D.  Not used if NW=1.
C
C  B       (input) DOUBLE PRECISION array, dimension (LDB,NW)
C          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
C          complex), column 1 contains the real part of B and column 2
C          contains the imaginary part.
C
C  LDB     (input) INTEGER
C          The leading dimension of B.  It must be at least NA.
C
C  WR      (input) DOUBLE PRECISION
C          The real part of the scalar "w".
C
C  WI      (input) DOUBLE PRECISION
C          The imaginary part of the scalar "w".  Not used if NW=1.
C
C  X       (output) DOUBLE PRECISION array, dimension (LDX,NW)
C          The NA x NW matrix X (unknowns), as computed by DLALN2.
C          If NW=2 ("w" is complex), on exit, column 1 will contain
C          the real part of X and column 2 will contain the imaginary
C          part.
C
C  LDX     (input) INTEGER
C          The leading dimension of X.  It must be at least NA.
C
C  SCALE   (output) DOUBLE PRECISION
C          The scale factor that B must be multiplied by to insure
C          that overflow does not occur when computing X.  Thus,
C          (ca A - w D) X  will be SCALE*B, not B (ignoring
C          perturbations of A.)  It will be at most 1.
C
C  XNORM   (output) DOUBLE PRECISION
C          The infinity-norm of X, when X is regarded as an NA x NW
C          real matrix.
C
C  INFO    (output) INTEGER
C          An error flag.  It will be set to zero if no error occurs,
C          a negative number if an argument is in error, or a positive
C          number if  ca A - w D  had to be perturbed.
C          The possible values are:
C          = 0: No error occurred, and (ca A - w D) did not have to be
C                 perturbed.
C          = 1: (ca A - w D) had to be perturbed to make its smallest
C               (or only) singular value greater than SMIN.
C          NOTE: In the interests of speed, this routine does not
C                check the inputs for errors.
C
C-----------------------------------------------------------------------
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
      INTEGER           INFO, LDA, LDB, LDX, NA, NW
      LOGICAL           LTRANS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), X(LDX,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21,
     *                  CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21,
     *                  LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R,
     *                  UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S,
     *                  UR22, XI1, XI2, XR1, XR2
      INTEGER           ICMAX, J
C     .. Local Arrays ..
      DOUBLE PRECISION  CI(2,2), CIV(4), CR(2,2), CRV(4)
      INTEGER           IPIVOT(4,4)
      LOGICAL           RSWAP(4), ZSWAP(4)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          A02ACF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Equivalences ..
      EQUIVALENCE       (CI(1,1),CIV(1)), (CR(1,1),CRV(1))
C     .. Data statements ..
      DATA              ZSWAP/.FALSE., .FALSE., .TRUE., .TRUE./
      DATA              RSWAP/.FALSE., .TRUE., .FALSE., .TRUE./
      DATA              IPIVOT/1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3,
     *                  2, 1/
C     .. Executable Statements ..
C
C     Compute BIGNUM
C
      SMLNUM = TWO*X02AMF()
      BIGNUM = ONE/SMLNUM
      SMINI = MAX(SMIN,SMLNUM)
C
C     Don't check for input errors
C
      INFO = 0
C
C     Standard Initializations
C
      SCALE = ONE
C
      IF (NA.EQ.1) THEN
C
C        1 x 1  (i.e., scalar) system   C X = B
C
         IF (NW.EQ.1) THEN
C
C           Real 1x1 system.
C
C           C = ca A - w D
C
            CSR = CA*A(1,1) - WR*D1
            CNORM = ABS(CSR)
C
C           If | C | < SMINI, use C = SMINI
C
            IF (CNORM.LT.SMINI) THEN
               CSR = SMINI
               CNORM = SMINI
               INFO = 1
            END IF
C
C           Check scaling for  X = B / C
C
            BNORM = ABS(B(1,1))
            IF (CNORM.LT.ONE .AND. BNORM.GT.ONE) THEN
               IF (BNORM.GT.BIGNUM*CNORM) SCALE = ONE/BNORM
            END IF
C
C           Compute X
C
            X(1,1) = (B(1,1)*SCALE)/CSR
            XNORM = ABS(X(1,1))
         ELSE
C
C           Complex 1x1 system (w is complex)
C
C           C = ca A - w D
C
            CSR = CA*A(1,1) - WR*D1
            CSI = -WI*D1
            CNORM = ABS(CSR) + ABS(CSI)
C
C           If | C | < SMINI, use C = SMINI
C
            IF (CNORM.LT.SMINI) THEN
               CSR = SMINI
               CSI = ZERO
               CNORM = SMINI
               INFO = 1
            END IF
C
C           Check scaling for  X = B / C
C
            BNORM = ABS(B(1,1)) + ABS(B(1,2))
            IF (CNORM.LT.ONE .AND. BNORM.GT.ONE) THEN
               IF (BNORM.GT.BIGNUM*CNORM) SCALE = ONE/BNORM
            END IF
C
C           Compute X
C
            CALL A02ACF(SCALE*B(1,1),SCALE*B(1,2),CSR,CSI,X(1,1),X(1,2))
            XNORM = ABS(X(1,1)) + ABS(X(1,2))
         END IF
C
      ELSE
C
C        2x2 System
C
C        Compute the real part of  C = ca A - w D  (or  ca A' - w D )
C
         CR(1,1) = CA*A(1,1) - WR*D1
         CR(2,2) = CA*A(2,2) - WR*D2
         IF (LTRANS) THEN
            CR(1,2) = CA*A(2,1)
            CR(2,1) = CA*A(1,2)
         ELSE
            CR(2,1) = CA*A(2,1)
            CR(1,2) = CA*A(1,2)
         END IF
C
         IF (NW.EQ.1) THEN
C
C           Real 2x2 system  (w is real)
C
C           Find the largest entry in C
C
            CMAX = ZERO
            ICMAX = 0
C
            DO 20 J = 1, 4
               IF (ABS(CRV(J)).GT.CMAX) THEN
                  CMAX = ABS(CRV(J))
                  ICMAX = J
               END IF
   20       CONTINUE
C
C           If norm(C) < SMINI, use SMINI*identity.
C
            IF (CMAX.LT.SMINI) THEN
               BNORM = MAX(ABS(B(1,1)),ABS(B(2,1)))
               IF (SMINI.LT.ONE .AND. BNORM.GT.ONE) THEN
                  IF (BNORM.GT.BIGNUM*SMINI) SCALE = ONE/BNORM
               END IF
               TEMP = SCALE/SMINI
               X(1,1) = TEMP*B(1,1)
               X(2,1) = TEMP*B(2,1)
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
C
C           Gaussian elimination with complete pivoting.
C
            UR11 = CRV(ICMAX)
            CR21 = CRV(IPIVOT(2,ICMAX))
            UR12 = CRV(IPIVOT(3,ICMAX))
            CR22 = CRV(IPIVOT(4,ICMAX))
            UR11R = ONE/UR11
            LR21 = UR11R*CR21
            UR22 = CR22 - UR12*LR21
C
C           If smaller pivot < SMINI, use SMINI
C
            IF (ABS(UR22).LT.SMINI) THEN
               UR22 = SMINI
               INFO = 1
            END IF
            IF (RSWAP(ICMAX)) THEN
               BR1 = B(2,1)
               BR2 = B(1,1)
            ELSE
               BR1 = B(1,1)
               BR2 = B(2,1)
            END IF
            BR2 = BR2 - LR21*BR1
            BBND = MAX(ABS(BR1*(UR22*UR11R)),ABS(BR2))
            IF (BBND.GT.ONE .AND. ABS(UR22).LT.ONE) THEN
               IF (BBND.GE.BIGNUM*ABS(UR22)) SCALE = ONE/BBND
            END IF
C
            XR2 = (BR2*SCALE)/UR22
            XR1 = (SCALE*BR1)*UR11R - XR2*(UR11R*UR12)
            IF (ZSWAP(ICMAX)) THEN
               X(1,1) = XR2
               X(2,1) = XR1
            ELSE
               X(1,1) = XR1
               X(2,1) = XR2
            END IF
            XNORM = MAX(ABS(XR1),ABS(XR2))
C
C           Further scaling if  norm(A) norm(X) > overflow
C
            IF (XNORM.GT.ONE .AND. CMAX.GT.ONE) THEN
               IF (XNORM.GT.BIGNUM/CMAX) THEN
                  TEMP = CMAX/BIGNUM
                  X(1,1) = TEMP*X(1,1)
                  X(2,1) = TEMP*X(2,1)
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         ELSE
C
C           Complex 2x2 system  (w is complex)
C
C           Find the largest entry in C
C
            CI(1,1) = -WI*D1
            CI(2,1) = ZERO
            CI(1,2) = ZERO
            CI(2,2) = -WI*D2
            CMAX = ZERO
            ICMAX = 0
C
            DO 40 J = 1, 4
               IF (ABS(CRV(J))+ABS(CIV(J)).GT.CMAX) THEN
                  CMAX = ABS(CRV(J)) + ABS(CIV(J))
                  ICMAX = J
               END IF
   40       CONTINUE
C
C           If norm(C) < SMINI, use SMINI*identity.
C
            IF (CMAX.LT.SMINI) THEN
               BNORM = MAX(ABS(B(1,1))+ABS(B(1,2)),ABS(B(2,1))
     *                 +ABS(B(2,2)))
               IF (SMINI.LT.ONE .AND. BNORM.GT.ONE) THEN
                  IF (BNORM.GT.BIGNUM*SMINI) SCALE = ONE/BNORM
               END IF
               TEMP = SCALE/SMINI
               X(1,1) = TEMP*B(1,1)
               X(2,1) = TEMP*B(2,1)
               X(1,2) = TEMP*B(1,2)
               X(2,2) = TEMP*B(2,2)
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
C
C           Gaussian elimination with complete pivoting.
C
            UR11 = CRV(ICMAX)
            UI11 = CIV(ICMAX)
            CR21 = CRV(IPIVOT(2,ICMAX))
            CI21 = CIV(IPIVOT(2,ICMAX))
            UR12 = CRV(IPIVOT(3,ICMAX))
            UI12 = CIV(IPIVOT(3,ICMAX))
            CR22 = CRV(IPIVOT(4,ICMAX))
            CI22 = CIV(IPIVOT(4,ICMAX))
            IF (ICMAX.EQ.1 .OR. ICMAX.EQ.4) THEN
C
C              Code when off-diagonals of pivoted C are real
C
               IF (ABS(UR11).GT.ABS(UI11)) THEN
                  TEMP = UI11/UR11
                  UR11R = ONE/(UR11*(ONE+TEMP**2))
                  UI11R = -TEMP*UR11R
               ELSE
                  TEMP = UR11/UI11
                  UI11R = -ONE/(UI11*(ONE+TEMP**2))
                  UR11R = -TEMP*UI11R
               END IF
               LR21 = CR21*UR11R
               LI21 = CR21*UI11R
               UR12S = UR12*UR11R
               UI12S = UR12*UI11R
               UR22 = CR22 - UR12*LR21
               UI22 = CI22 - UR12*LI21
            ELSE
C
C              Code when diagonals of pivoted C are real
C
               UR11R = ONE/UR11
               UI11R = ZERO
               LR21 = CR21*UR11R
               LI21 = CI21*UR11R
               UR12S = UR12*UR11R
               UI12S = UI12*UR11R
               UR22 = CR22 - UR12*LR21 + UI12*LI21
               UI22 = -UR12*LI21 - UI12*LR21
            END IF
            U22ABS = ABS(UR22) + ABS(UI22)
C
C           If smaller pivot < SMINI, use SMINI
C
            IF (U22ABS.LT.SMINI) THEN
               UR22 = SMINI
               UI22 = ZERO
               INFO = 1
            END IF
            IF (RSWAP(ICMAX)) THEN
               BR2 = B(1,1)
               BR1 = B(2,1)
               BI2 = B(1,2)
               BI1 = B(2,2)
            ELSE
               BR1 = B(1,1)
               BR2 = B(2,1)
               BI1 = B(1,2)
               BI2 = B(2,2)
            END IF
            BR2 = BR2 - LR21*BR1 + LI21*BI1
            BI2 = BI2 - LI21*BR1 - LR21*BI1
            BBND = MAX((ABS(BR1)+ABS(BI1))*(U22ABS*(ABS(UR11R)
     *             +ABS(UI11R))),ABS(BR2)+ABS(BI2))
            IF (BBND.GT.ONE .AND. U22ABS.LT.ONE) THEN
               IF (BBND.GE.BIGNUM*U22ABS) THEN
                  SCALE = ONE/BBND
                  BR1 = SCALE*BR1
                  BI1 = SCALE*BI1
                  BR2 = SCALE*BR2
                  BI2 = SCALE*BI2
               END IF
            END IF
C
            CALL A02ACF(BR2,BI2,UR22,UI22,XR2,XI2)
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2
            IF (ZSWAP(ICMAX)) THEN
               X(1,1) = XR2
               X(2,1) = XR1
               X(1,2) = XI2
               X(2,2) = XI1
            ELSE
               X(1,1) = XR1
               X(2,1) = XR2
               X(1,2) = XI1
               X(2,2) = XI2
            END IF
            XNORM = MAX(ABS(XR1)+ABS(XI1),ABS(XR2)+ABS(XI2))
C
C           Further scaling if  norm(A) norm(X) > overflow
C
            IF (XNORM.GT.ONE .AND. CMAX.GT.ONE) THEN
               IF (XNORM.GT.BIGNUM/CMAX) THEN
                  TEMP = CMAX/BIGNUM
                  X(1,1) = TEMP*X(1,1)
                  X(2,1) = TEMP*X(2,1)
                  X(1,2) = TEMP*X(1,2)
                  X(2,2) = TEMP*X(2,2)
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F08QHX (DLALN2)
C
      END
