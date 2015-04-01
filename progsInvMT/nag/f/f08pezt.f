      SUBROUTINE F08PEZ(WANTT,WANTZ,N,ILO,IHI,H,LDH,WR,WI,ILOZ,IHIZ,Z,
     *                  LDZ,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAHQR(WANTT,WANTZ,N,ILO,IHI,H,LDH,WR,WI,ILOZ,
C    *                  IHIZ,Z,LDZ,INFO)
C
C  Purpose
C  =======
C
C  DLAHQR is an auxiliary routine called by DHSEQR to update the
C  eigenvalues and Schur decomposition already computed by DHSEQR, by
C  dealing with the Hessenberg submatrix in rows and columns ILO to IHI.
C
C  Arguments
C  =========
C
C  WANTT   (input) LOGICAL
C          = .TRUE. : the full Schur form T is required;
C          = .FALSE.: only eigenvalues are required.
C
C  WANTZ   (input) LOGICAL
C          = .TRUE. : the matrix of Schur vectors Z is required;
C          = .FALSE.: Schur vectors are not required.
C
C  N       (input) INTEGER
C          The order of the matrix H.  N >= 0.
C
C  ILO     (input) INTEGER
C  IHI     (input) INTEGER
C          It is assumed that H is already upper quasi-triangular in
C          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
C          ILO = 1). DLAHQR works primarily with the Hessenberg
C          submatrix in rows and columns ILO to IHI, but applies
C          transformations to all of H if WANTT is .TRUE..
C          1 <= ILO <= max(1,IHI); IHI <= N.
C
C  H       (input/output) REAL array, dimension (LDH,N)
C          On entry, the upper Hessenberg matrix H.
C          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
C          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
C          standard form. If WANTT is .FALSE., the contents of H are
C          unspecified on exit.
C
C  LDH     (input) INTEGER
C          The leading dimension of the array H. LDH >= max(1,N).
C
C  WR      (output) REAL array, dimension (N)
C  WI      (output) REAL array, dimension (N)
C          The real and imaginary parts, respectively, of the computed
C          eigenvalues ILO to IHI are stored in the corresponding
C          elements of WR and WI. If two eigenvalues are computed as a
C          complex conjugate pair, they are stored in consecutive
C          elements of WR and WI, say the i-th and (i+1)th, with
C          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
C          eigenvalues are stored in the same order as on the diagonal
C          of the Schur form returned in H, with WR(i) = H(i,i), and, if
C          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
C          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
C
C  ILOZ    (input) INTEGER
C  IHIZ    (input) INTEGER
C          Specify the rows of Z to which transformations must be
C          applied if WANTZ is .TRUE..
C          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
C
C  Z       (input/output) REAL array, dimension (LDZ,N)
C          If WANTZ is .TRUE., on entry Z must contain the current
C          matrix Z of transformations accumulated by DHSEQR, and on
C          exit Z has been updated; transformations are applied only to
C          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
C          If WANTZ is .FALSE., Z is not referenced.
C
C  LDZ     (input) INTEGER
C          The leading dimension of the array Z. LDZ >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          > 0: DLAHQR failed to compute all the eigenvalues ILO to IHI
C               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
C               elements i+1:ihi of WR and WI contain those eigenvalues
C               which have been successfully computed.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  DAT1, DAT2
      PARAMETER         (DAT1=0.75D+0,DAT2=-0.4375D+0)
C     .. Scalar Arguments ..
      INTEGER           IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL           WANTT, WANTZ
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LDH,*), WI(*), WR(*), Z(LDZ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, H00, H10, H11, H12, H21, H22, H33, H33S,
     *                  H43H34, H44, H44S, OVFL, S, SMLNUM, SN, SUM, T1,
     *                  T2, T3, TST1, ULP, UNFL, V1, V2, V3
      INTEGER           I, I1, I2, ITN, ITS, J, K, L, M, NH, NR, NZ
C     .. Local Arrays ..
      DOUBLE PRECISION  V(3), WORK(1)
C     .. External Functions ..
      DOUBLE PRECISION  F06RMF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          F06RMF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DROT, F08AEV, F08PEY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
      IF (ILO.EQ.IHI) THEN
         WR(ILO) = H(ILO,ILO)
         WI(ILO) = ZERO
         RETURN
      END IF
C
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
C
C     Set machine-dependent constants for the stopping criterion.
C     If norm(H) <= sqrt(OVFL), overflow should not occur.
C
      UNFL = X02AMF()
      OVFL = ONE/UNFL
      ULP = X02AJF()*X02BHF()
      SMLNUM = UNFL*(NH/ULP)
C
C     I1 and I2 are the indices of the first row and last column of H
C     to which transformations must be applied. If eigenvalues only are
C     being computed, I1 and I2 are set inside the main loop.
C
      IF (WANTT) THEN
         I1 = 1
         I2 = N
      END IF
C
C     ITN is the total number of QR iterations allowed.
C
      ITN = 30*NH
C
C     The main loop begins here. I is the loop index and decreases from
C     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
C     with the active submatrix in rows and columns L to I.
C     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
C     H(L,L-1) is negligible so that the matrix splits.
C
      I = IHI
   20 CONTINUE
      L = ILO
      IF (I.LT.ILO) GO TO 300
C
C     Perform QR iterations on rows and columns ILO to I until a
C     submatrix of order 1 or 2 splits off at the bottom because a
C     subdiagonal element has become negligible.
C
      DO 260 ITS = 0, ITN
C
C        Look for a single small subdiagonal element.
C
         DO 40 K = I, L + 1, -1
            TST1 = ABS(H(K-1,K-1)) + ABS(H(K,K))
            IF (TST1.EQ.ZERO) TST1 = F06RMF('1',I-L+1,H(L,L),LDH,WORK)
            IF (ABS(H(K,K-1)).LE.MAX(ULP*TST1,SMLNUM)) GO TO 60
   40    CONTINUE
   60    CONTINUE
         L = K
         IF (L.GT.ILO) THEN
C
C           H(L,L-1) is negligible
C
            H(L,L-1) = ZERO
         END IF
C
C        Exit from loop if a submatrix of order 1 or 2 has split off.
C
         IF (L.GE.I-1) GO TO 280
C
C        Now the active submatrix is in rows and columns L to I. If
C        eigenvalues only are being computed, only the active submatrix
C        need be transformed.
C
         IF ( .NOT. WANTT) THEN
            I1 = L
            I2 = I
         END IF
C
         IF (ITS.EQ.10 .OR. ITS.EQ.20) THEN
C
C           Exceptional shift.
C
            S = ABS(H(I,I-1)) + ABS(H(I-1,I-2))
            H44 = DAT1*S
            H33 = H44
            H43H34 = DAT2*S*S
         ELSE
C
C           Prepare to use Wilkinson's double shift
C
            H44 = H(I,I)
            H33 = H(I-1,I-1)
            H43H34 = H(I,I-1)*H(I-1,I)
         END IF
C
C        Look for two consecutive small subdiagonal elements.
C
         DO 80 M = I - 2, L, -1
C
C           Determine the effect of starting the double-shift QR
C           iteration at row M, and see if this would make H(M,M-1)
C           negligible.
C
            H11 = H(M,M)
            H22 = H(M+1,M+1)
            H21 = H(M+1,M)
            H12 = H(M,M+1)
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = (H33S*H44S-H43H34)/H21 + H12
            V2 = H22 - H11 - H33S - H44S
            V3 = H(M+2,M+1)
            S = ABS(V1) + ABS(V2) + ABS(V3)
            V1 = V1/S
            V2 = V2/S
            V3 = V3/S
            V(1) = V1
            V(2) = V2
            V(3) = V3
            IF (M.EQ.L) GO TO 100
            H00 = H(M-1,M-1)
            H10 = H(M,M-1)
            TST1 = ABS(V1)*(ABS(H00)+ABS(H11)+ABS(H22))
            IF (ABS(H10)*(ABS(V2)+ABS(V3)).LE.ULP*TST1) GO TO 100
   80    CONTINUE
  100    CONTINUE
C
C        Double-shift QR step
C
         DO 240 K = M, I - 1
C
C           The first iteration of this loop determines a reflection G
C           from the vector V and applies it from left and right to H,
C           thus creating a nonzero bulge below the subdiagonal.
C
C           Each subsequent iteration determines a reflection G to
C           restore the Hessenberg form in the (K-1)th column, and thus
C           chases the bulge one step toward the bottom of the active
C           submatrix. NR is the order of G.
C
            NR = MIN(3,I-K+1)
            IF (K.GT.M) CALL DCOPY(NR,H(K,K-1),1,V,1)
            CALL F08AEV(NR,V(1),V(2),1,T1)
            IF (K.GT.M) THEN
               H(K,K-1) = V(1)
               H(K+1,K-1) = ZERO
               IF (K.LT.I-1) H(K+2,K-1) = ZERO
            ELSE IF (M.GT.L) THEN
               H(K,K-1) = -H(K,K-1)
            END IF
            V2 = V(2)
            T2 = T1*V2
            IF (NR.EQ.3) THEN
               V3 = V(3)
               T3 = T1*V3
C
C              Apply G from the left to transform the rows of the matrix
C              in columns K to I2.
C
               DO 120 J = K, I2
                  SUM = H(K,J) + V2*H(K+1,J) + V3*H(K+2,J)
                  H(K,J) = H(K,J) - SUM*T1
                  H(K+1,J) = H(K+1,J) - SUM*T2
                  H(K+2,J) = H(K+2,J) - SUM*T3
  120          CONTINUE
C
C              Apply G from the right to transform the columns of the
C              matrix in rows I1 to min(K+3,I).
C
               DO 140 J = I1, MIN(K+3,I)
                  SUM = H(J,K) + V2*H(J,K+1) + V3*H(J,K+2)
                  H(J,K) = H(J,K) - SUM*T1
                  H(J,K+1) = H(J,K+1) - SUM*T2
                  H(J,K+2) = H(J,K+2) - SUM*T3
  140          CONTINUE
C
               IF (WANTZ) THEN
C
C                 Accumulate transformations in the matrix Z
C
                  DO 160 J = ILOZ, IHIZ
                     SUM = Z(J,K) + V2*Z(J,K+1) + V3*Z(J,K+2)
                     Z(J,K) = Z(J,K) - SUM*T1
                     Z(J,K+1) = Z(J,K+1) - SUM*T2
                     Z(J,K+2) = Z(J,K+2) - SUM*T3
  160             CONTINUE
               END IF
            ELSE IF (NR.EQ.2) THEN
C
C              Apply G from the left to transform the rows of the matrix
C              in columns K to I2.
C
               DO 180 J = K, I2
                  SUM = H(K,J) + V2*H(K+1,J)
                  H(K,J) = H(K,J) - SUM*T1
                  H(K+1,J) = H(K+1,J) - SUM*T2
  180          CONTINUE
C
C              Apply G from the right to transform the columns of the
C              matrix in rows I1 to min(K+3,I).
C
               DO 200 J = I1, I
                  SUM = H(J,K) + V2*H(J,K+1)
                  H(J,K) = H(J,K) - SUM*T1
                  H(J,K+1) = H(J,K+1) - SUM*T2
  200          CONTINUE
C
               IF (WANTZ) THEN
C
C                 Accumulate transformations in the matrix Z
C
                  DO 220 J = ILOZ, IHIZ
                     SUM = Z(J,K) + V2*Z(J,K+1)
                     Z(J,K) = Z(J,K) - SUM*T1
                     Z(J,K+1) = Z(J,K+1) - SUM*T2
  220             CONTINUE
               END IF
            END IF
  240    CONTINUE
C
  260 CONTINUE
C
C     Failure to converge in remaining number of iterations
C
      INFO = I
      RETURN
C
  280 CONTINUE
C
      IF (L.EQ.I) THEN
C
C        H(I,I-1) is negligible: one eigenvalue has converged.
C
         WR(I) = H(I,I)
         WI(I) = ZERO
      ELSE IF (L.EQ.I-1) THEN
C
C        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
C
C        Transform the 2-by-2 submatrix to standard Schur form,
C        and compute and store the eigenvalues.
C
         CALL F08PEY(H(I-1,I-1),H(I-1,I),H(I,I-1),H(I,I),WR(I-1),WI(I-1)
     *               ,WR(I),WI(I),CS,SN)
C
         IF (WANTT) THEN
C
C           Apply the transformation to the rest of H.
C
            IF (I2.GT.I) CALL DROT(I2-I,H(I-1,I+1),LDH,H(I,I+1),LDH,CS,
     *                             SN)
            CALL DROT(I-I1-1,H(I1,I-1),1,H(I1,I),1,CS,SN)
         END IF
         IF (WANTZ) THEN
C
C           Apply the transformation to Z.
C
            CALL DROT(NZ,Z(ILOZ,I-1),1,Z(ILOZ,I),1,CS,SN)
         END IF
      END IF
C
C     Decrement number of remaining iterations, and return to start of
C     the main loop with new value of I.
C
      ITN = ITN - ITS
      I = L - 1
      GO TO 20
C
  300 CONTINUE
      RETURN
C
C     End of F08PEZ (DLAHQR)
C
      END
