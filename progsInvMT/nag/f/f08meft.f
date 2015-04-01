      SUBROUTINE F08MEF(UPLO,N,NCVT,NRU,NCC,D,E,VT,LDVT,U,LDU,C,LDC,
     *                  WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 RE-ISSUE. NAG COPYRIGHT 1995.
C     .. Entry Points ..
      ENTRY             DBDSQR(UPLO,N,NCVT,NRU,NCC,D,E,VT,LDVT,U,LDU,C,
     *                  LDC,WORK,INFO)
C
C  Purpose
C  =======
C
C  DBDSQR computes the singular value decomposition (SVD) of
C  a real N-by-N (upper or lower) bidiagonal matrix B:
C  B = Q * S * P' (P' denotes the transpose of P), where S is a
C  diagonal matrix with non-negative diagonal elements (the singular
C  values of B), and Q and P are orthogonal matrices.
C
C  The routine computes S, and optionally computes U * Q, P' * VT,
C  or Q' * C, for given real input matrices U, VT, and C.
C
C  See "Computing  Small Singular Values of Bidiagonal Matrices With
C  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
C  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
C  no. 5, pp. 873-912, Sept 1990),
C
C  "Accurate singular values and differential qd algorithms," by
C  K. V. Fernando and B. N. Parlett,
C  Numer. Math., Vol-67, No. 2, pp. 191-230,1994.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          = 'U':  B is upper bidiagonal;
C          = 'L':  B is lower bidiagonal.
C
C  N       (input) INTEGER
C          The order of the matrix B.  N >= 0.
C
C  NCVT    (input) INTEGER
C          The number of columns of the matrix VT. NCVT >= 0.
C
C  NRU     (input) INTEGER
C          The number of rows of the matrix U. NRU >= 0.
C
C  NCC     (input) INTEGER
C          The number of columns of the matrix C. NCC >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, the n diagonal elements of the bidiagonal matrix B.
C          On exit, if INFO=0, the singular values of B in decreasing
C          order.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          On entry, the elements of E contain the
C          offdiagonal elements of the bidiagonal matrix whose SVD
C          is desired. On normal exit (INFO = 0), E is destroyed.
C          If the algorithm does not converge (INFO > 0), D and E
C          will contain the diagonal and superdiagonal elements of a
C          bidiagonal matrix orthogonally equivalent to the one given
C          as input.
C
C  VT      (input/output) DOUBLE PRECISION array, dimension (LDVT,NCVT)
C          On entry, an N-by-NCVT matrix VT.
C          On exit, VT is overwritten by P' * VT.
C          VT is not referenced if NCVT = 0.
C
C  LDVT    (input) INTEGER
C          The leading dimension of the array VT.
C          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
C
C  U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C          On entry, an NRU-by-N matrix U.
C          On exit, U is overwritten by U * Q.
C          U is not referenced if NRU = 0.
C
C  LDU     (input) INTEGER
C          The leading dimension of the array U.  LDU >= max(1,NRU).
C
C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,NCC)
C          On entry, an N-by-NCC matrix C.
C          On exit, C is overwritten by Q' * C.
C          C is not referenced if NCC = 0.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C.
C          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N) if only
C          singular values are wanted (NCVT = NRU = NCC = 0), and
C          dimension (max(1,4*N-4)) otherwise.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  If INFO = -i, the i-th argument had an illegal value
C          > 0:  the algorithm did not converge; D and E contain the
C                elements of a bidiagonal matrix which is orthogonally
C                similar to the input matrix B;  if INFO = i, i
C                elements of E have not converged to zero.
C
C  Internal Parameters
C  ===================
C
C  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
C          TOLMUL controls the convergence criterion of the QR loop.
C          If it is positive, TOLMUL*EPS is the desired relative
C             precision in the computed singular values.
C          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
C             desired absolute accuracy in the computed singular
C             values (corresponds to relative accuracy
C             abs(TOLMUL*EPS) in the largest singular value).
C          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
C             between 10 (for fast convergence) and .1/EPS
C             (for there to be some accuracy in the results).
C          Default is to lose at most either one eighth or 2 of the
C             available decimal digits in each computed singular value
C             (whichever is smaller).
C
C  MAXITR  INTEGER, default = 6
C          MAXITR controls the maximum number of passes of the
C          algorithm through its inner loop. The algorithms stops
C          (and so fails to converge) if the number of passes
C          through the inner loop exceeds MAXITR*N**2.
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      DOUBLE PRECISION  NEGONE
      PARAMETER         (NEGONE=-1.0D0)
      DOUBLE PRECISION  HNDRTH
      PARAMETER         (HNDRTH=0.01D0)
      DOUBLE PRECISION  TEN
      PARAMETER         (TEN=10.0D0)
      DOUBLE PRECISION  HNDRD
      PARAMETER         (HNDRD=100.0D0)
      DOUBLE PRECISION  MEIGTH
      PARAMETER         (MEIGTH=-0.125D0)
      INTEGER           MAXITR
      PARAMETER         (MAXITR=6)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,*), D(*), E(*), U(LDU,*), VT(LDVT,*),
     *                  WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSE, ABSS, COSL, COSR, EPS, F, G, H, MU, R,
     *                  SHIFT, SIGMN, SIGMX, SINL, SINR, SLL, SMAX,
     *                  SMIN, SMINL, SMINOA, THRESH, TOL, TOLMUL, UNFL
      INTEGER           I, ISUB, ITER, J, LAST, LL, LLL, M, MAXIT, NM1,
     *                  NM12, NM13, OLDLL, OLDM
      LOGICAL           FORWRD, UPPER
      CHARACTER         DIRECT
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DROT, DSCAL, DSWAP, F06AAZ, F06QXF, F08HEW,
     *                  F08MEW, F08MEX, F08MEY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SIGN, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = UPLO .EQ. 'U' .OR. UPLO .EQ. 'u'
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (NCVT.LT.0) THEN
         INFO = -3
      ELSE IF (NRU.LT.0) THEN
         INFO = -4
      ELSE IF (NCC.LT.0) THEN
         INFO = -5
      ELSE IF ((NCVT.EQ.0 .AND. LDVT.LT.1)
     *         .OR. (NCVT.GT.0 .AND. LDVT.LT.MAX(1,N))) THEN
         INFO = -9
      ELSE IF (LDU.LT.MAX(1,NRU)) THEN
         INFO = -11
      ELSE IF ((NCC.EQ.0 .AND. LDC.LT.1)
     *         .OR. (NCC.GT.0 .AND. LDC.LT.MAX(1,N))) THEN
         INFO = -13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08MEF/DBDSQR',-INFO)
         RETURN
      END IF
      IF (N.EQ.0) RETURN
      IF (N.EQ.1) GO TO 300
C
C     If no singular vectors desired, use qd algorithm
C
      IF (NCVT.EQ.0 .AND. NRU.EQ.0 .AND. NCC.EQ.0) THEN
         CALL F08MEW(N,D,E,WORK,INFO)
         RETURN
      END IF
C
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
C
C     Get machine constants
C
      EPS = X02AJF()
      UNFL = X02AMF()
      TOLMUL = MAX(TEN,MIN(HNDRD,EPS**MEIGTH))
      TOL = TOLMUL*EPS
C
C     If matrix lower bidiagonal, rotate to be upper bidiagonal
C     by applying Givens rotations on the left
C
      IF ( .NOT. UPPER) THEN
         DO 20 I = 1, N - 1
            CALL F08HEW(D(I),E(I),COSR,SINR,R)
            D(I) = R
            E(I) = SINR*D(I+1)
            D(I+1) = COSR*D(I+1)
            WORK(I) = COSR
            WORK(NM1+I) = SINR
   20    CONTINUE
C
C        Update singular vectors if desired
C
         IF (NRU.GT.0) CALL F06QXF('R','V','F',NRU,N,1,N,WORK(1),WORK(N)
     *                             ,U,LDU)
         IF (NCC.GT.0) CALL F06QXF('L','V','F',N,NCC,1,N,WORK(1),WORK(N)
     *                             ,C,LDC)
      END IF
C
C     Compute approximate maximum, minimum singular values
C
      SMAX = ABS(D(N))
      DO 40 I = 1, N - 1
         SMAX = MAX(SMAX,ABS(D(I)),ABS(E(I)))
   40 CONTINUE
      SMINL = ZERO
      IF (TOL.GE.ZERO) THEN
         SMINOA = ABS(D(1))
         IF (SMINOA.EQ.ZERO) GO TO 80
         MU = SMINOA
         DO 60 I = 2, N
            MU = ABS(D(I))*(MU/(MU+ABS(E(I-1))))
            SMINOA = MIN(SMINOA,MU)
            IF (SMINOA.EQ.ZERO) GO TO 80
   60    CONTINUE
   80    CONTINUE
         SMINOA = SMINOA/SQRT(DBLE(N))
      END IF
C
C     Prepare for main iteration loop for the singular values
C
      MAXIT = MAXITR*N*N
      ITER = 0
      OLDLL = -1
      OLDM = -1
      IF (TOL.GE.ZERO) THEN
C
C        Relative accuracy desired
C
         THRESH = MAX(TOL*SMINOA,MAXIT*UNFL)
      ELSE
C
C        Absolute accuracy desired
C
         THRESH = MAX(ABS(TOL)*SMAX,MAXIT*UNFL)
      END IF
C
C     M points to last entry of unconverged part of matrix
C
      M = N
C
C     Begin main iteration loop
C
  100 CONTINUE
C
C     Check for convergence or exceeding iteration count
C
      IF (M.LE.1) GO TO 300
      IF (ITER.GT.MAXIT) GO TO 380
C
C     Find diagonal block of matrix to work on
C
      IF (TOL.LT.ZERO .AND. ABS(D(M)).LE.THRESH) D(M) = ZERO
      SMAX = ABS(D(M))
      SMIN = SMAX
      DO 120 LLL = 1, M
         LL = M - LLL
         IF (LL.EQ.0) GO TO 160
         ABSS = ABS(D(LL))
         ABSE = ABS(E(LL))
         IF (TOL.LT.ZERO .AND. ABSS.LE.THRESH) D(LL) = ZERO
         IF (ABSE.LE.THRESH) GO TO 140
         SMIN = MIN(SMIN,ABSS)
         SMAX = MAX(SMAX,ABSS,ABSE)
  120 CONTINUE
  140 CONTINUE
      E(LL) = ZERO
C
C     Matrix splits since E(LL) = 0
C
      IF (LL.EQ.M-1) THEN
C
C        Convergence of bottom singular value, return to top of loop
C
         M = M - 1
         GO TO 100
      END IF
  160 CONTINUE
      LL = LL + 1
C
C     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
C
      IF (LL.EQ.M-1) THEN
C
C        2 by 2 block, handle separately
C
         CALL F08MEY(D(M-1),E(M-1),D(M),SIGMN,SIGMX,SINR,COSR,SINL,COSL)
         D(M-1) = SIGMX
         E(M-1) = ZERO
         D(M) = SIGMN
C
C        Compute singular vectors, if desired
C
         IF (NCVT.GT.0) CALL DROT(NCVT,VT(M-1,1),LDVT,VT(M,1),LDVT,COSR,
     *                            SINR)
         IF (NRU.GT.0) CALL DROT(NRU,U(1,M-1),1,U(1,M),1,COSL,SINL)
         IF (NCC.GT.0) CALL DROT(NCC,C(M-1,1),LDC,C(M,1),LDC,COSL,SINL)
         M = M - 2
         GO TO 100
      END IF
C
C     If working on new submatrix, choose shift direction
C     (from larger end diagonal entry towards smaller)
C
      IF (LL.GT.OLDM .OR. M.LT.OLDLL) FORWRD = ABS(D(LL)) .GE. ABS(D(M))
C
C     Apply convergence tests
C
      IF (FORWRD) THEN
C
C        Run convergence test in forward direction
C        First apply standard test to bottom of matrix
C
         IF (ABS(E(M-1)).LE.ABS(TOL)*ABS(D(M))
     *       .OR. (TOL.LT.ZERO .AND. ABS(E(M-1)).LE.THRESH)) THEN
            E(M-1) = ZERO
            GO TO 100
         END IF
C
         IF (TOL.GE.ZERO) THEN
C
C           If relative accuracy desired,
C           apply convergence criterion forward
C
            MU = ABS(D(LL))
            SMINL = MU
            DO 180 LLL = LL, M - 1
               IF (ABS(E(LLL)).LE.TOL*MU) THEN
                  E(LLL) = ZERO
                  GO TO 100
               END IF
               MU = ABS(D(LLL+1))*(MU/(MU+ABS(E(LLL))))
               SMINL = MIN(SMINL,MU)
  180       CONTINUE
         END IF
      ELSE
C
C        Run convergence test in backward direction
C        First apply standard test to top of matrix
C
         IF (ABS(E(LL)).LE.ABS(TOL)*ABS(D(LL))
     *       .OR. (TOL.LT.ZERO .AND. ABS(E(LL)).LE.THRESH)) THEN
            E(LL) = ZERO
            GO TO 100
         END IF
C
         IF (TOL.GE.ZERO) THEN
C
C           If relative accuracy desired,
C           apply convergence criterion backward
C
            MU = ABS(D(M))
            SMINL = MU
            DO 200 LLL = M - 1, LL, -1
               IF (ABS(E(LLL)).LE.TOL*MU) THEN
                  E(LLL) = ZERO
                  GO TO 100
               END IF
               MU = ABS(D(LLL))*(MU/(MU+ABS(E(LLL))))
               SMINL = MIN(SMINL,MU)
  200       CONTINUE
         END IF
      END IF
      OLDLL = LL
      OLDM = M
C
C     Compute shift.  First, test if shifting would ruin relative
C     accuracy, and if so set the shift to zero.
C
      IF (TOL.GE.ZERO .AND. N*TOL*(SMINL/SMAX).LE.MAX(EPS,HNDRTH*TOL))
     *    THEN
C
C        Use a zero shift to avoid loss of relative accuracy
C
         SHIFT = ZERO
      ELSE
C
C        Compute the shift from 2-by-2 block at end of matrix
C
         IF (FORWRD) THEN
            SLL = ABS(D(LL))
            CALL F08MEX(D(M-1),E(M-1),D(M),SHIFT,R)
         ELSE
            SLL = ABS(D(M))
            CALL F08MEX(D(LL),E(LL),D(LL+1),SHIFT,R)
         END IF
C
C        Test if shift negligible, and if so set to zero
C
         IF (SLL.GT.ZERO) THEN
            IF ((SHIFT/SLL)**2.LT.EPS) SHIFT = ZERO
         END IF
      END IF
C
C     Increment iteration count
C
      ITER = ITER + M - LL
C
      IF (FORWRD) THEN
C
C        Chase bulge from top to bottom
C
         DIRECT = 'F'
         LAST = M - 1
C
         IF (SHIFT.EQ.ZERO) THEN
C
C           If SHIFT = 0, do simplified QR iteration
C
            COSR = ONE
            COSL = ONE
            DO 220 I = LL, M - 1
               CALL F08HEW(D(I)*COSR,E(I),COSR,SINR,R)
               IF (I.GT.LL) E(I-1) = SINL*R
               CALL F08HEW(COSL*R,D(I+1)*SINR,COSL,SINL,D(I))
               WORK(I) = COSL
               WORK(I+NM1) = SINL
               WORK(I+NM12) = COSR
               WORK(I+NM13) = SINR
  220       CONTINUE
            H = D(M)*COSR
            D(M) = H*COSL
            E(M-1) = H*SINL
C
         ELSE
C
C           Use nonzero shift
C
            F = (ABS(D(LL))-SHIFT)*(SIGN(ONE,D(LL))+SHIFT/D(LL))
            G = E(LL)
            DO 240 I = LL, M - 1
               CALL F08HEW(F,G,COSR,SINR,R)
               IF (I.GT.LL) E(I-1) = R
               F = COSR*D(I) + SINR*E(I)
               E(I) = COSR*E(I) - SINR*D(I)
               G = SINR*D(I+1)
               D(I+1) = COSR*D(I+1)
               CALL F08HEW(F,G,COSL,SINL,R)
               D(I) = R
               F = COSL*E(I) + SINL*D(I+1)
               D(I+1) = COSL*D(I+1) - SINL*E(I)
               IF (I.LT.M-1) THEN
                  G = SINL*E(I+1)
                  E(I+1) = COSL*E(I+1)
               END IF
               WORK(I) = COSL
               WORK(I+NM1) = SINL
               WORK(I+NM12) = COSR
               WORK(I+NM13) = SINR
  240       CONTINUE
            E(M-1) = F
C
         END IF
      ELSE
C
C        Chase bulge from bottom to top
C
         DIRECT = 'B'
         LAST = LL
C
         IF (SHIFT.EQ.ZERO) THEN
C
C           If SHIFT = 0, do simplified QR iteration
C
            COSR = ONE
            COSL = ONE
            DO 260 I = M, LL + 1, -1
               CALL F08HEW(D(I)*COSR,E(I-1),COSR,SINR,R)
               IF (I.LT.M) E(I) = SINL*R
               CALL F08HEW(COSL*R,D(I-1)*SINR,COSL,SINL,D(I))
               WORK(I-1) = COSR
               WORK(I-1+NM1) = -SINR
               WORK(I-1+NM12) = COSL
               WORK(I-1+NM13) = -SINL
  260       CONTINUE
            H = D(LL)*COSR
            D(LL) = H*COSL
            E(LL) = H*SINL
C
         ELSE
C
C           Nonzero shift
C
            F = (ABS(D(M))-SHIFT)*(SIGN(ONE,D(M))+SHIFT/D(M))
            G = E(M-1)
            DO 280 I = M, LL + 1, -1
               CALL F08HEW(F,G,COSR,SINR,R)
               IF (I.LT.M) E(I) = R
               F = COSR*D(I) + SINR*E(I-1)
               E(I-1) = COSR*E(I-1) - SINR*D(I)
               G = SINR*D(I-1)
               D(I-1) = COSR*D(I-1)
               CALL F08HEW(F,G,COSL,SINL,R)
               D(I) = R
               F = COSL*E(I-1) + SINL*D(I-1)
               D(I-1) = COSL*D(I-1) - SINL*E(I-1)
               IF (I.GT.LL+1) THEN
                  G = SINL*E(I-2)
                  E(I-2) = COSL*E(I-2)
               END IF
               WORK(I-1) = COSR
               WORK(I-1+NM1) = -SINR
               WORK(I-1+NM12) = COSL
               WORK(I-1+NM13) = -SINL
  280       CONTINUE
            E(LL) = F
C
         END IF
      END IF
C
C     Test convergence
C
      IF (ABS(E(LAST)).LE.THRESH) E(LAST) = ZERO
C
C     Update singular vectors if desired
C
      IF (NCVT.GT.0) CALL F06QXF('L','V',DIRECT,M-LL+1,NCVT,1,M-LL+1,
     *                           WORK(NM12+LL),WORK(NM13+LL),VT(LL,1),
     *                           LDVT)
      IF (NRU.GT.0) CALL F06QXF('R','V',DIRECT,NRU,M-LL+1,1,M-LL+1,
     *                          WORK(LL),WORK(NM1+LL),U(1,LL),LDU)
      IF (NCC.GT.0) CALL F06QXF('L','V',DIRECT,M-LL+1,NCC,1,M-LL+1,
     *                          WORK(LL),WORK(NM1+LL),C(LL,1),LDC)
C
C     QR iteration finished, go back and check convergence
C
      GO TO 100
C
C     All singular values converged, so make them positive
C
  300 CONTINUE
      DO 320 I = 1, N
         IF (D(I).LT.ZERO) THEN
            D(I) = -D(I)
C
C           Change sign of singular vectors, if desired
C
            IF (NCVT.GT.0) CALL DSCAL(NCVT,NEGONE,VT(I,1),LDVT)
         END IF
  320 CONTINUE
C
C     Sort the singular values into decreasing order (insertion sort on
C     singular values, but only one transposition per singular vector)
C
      DO 360 I = 1, N - 1
C
C        Scan for smallest D(I)
C
         ISUB = 1
         SMIN = D(1)
         DO 340 J = 2, N + 1 - I
            IF (D(J).LE.SMIN) THEN
               ISUB = J
               SMIN = D(J)
            END IF
  340    CONTINUE
         IF (ISUB.NE.N+1-I) THEN
C
C           Swap singular values and vectors
C
            D(ISUB) = D(N+1-I)
            D(N+1-I) = SMIN
            IF (NCVT.GT.0) CALL DSWAP(NCVT,VT(ISUB,1),LDVT,VT(N+1-I,1),
     *                                LDVT)
            IF (NRU.GT.0) CALL DSWAP(NRU,U(1,ISUB),1,U(1,N+1-I),1)
            IF (NCC.GT.0) CALL DSWAP(NCC,C(ISUB,1),LDC,C(N+1-I,1),LDC)
         END IF
  360 CONTINUE
      GO TO 420
C
C     Maximum number of iterations exceeded, failure to converge
C
  380 CONTINUE
      INFO = 0
      DO 400 I = 1, N - 1
         IF (E(I).NE.ZERO) INFO = INFO + 1
  400 CONTINUE
  420 CONTINUE
      RETURN
C
C     End of F08MEF (DBDSQR)
C
      END
