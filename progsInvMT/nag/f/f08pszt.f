      SUBROUTINE F08PSZ(WANTT,WANTZ,N,ILO,IHI,H,LDH,W,ILOZ,IHIZ,Z,LDZ,
     *                  INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLAHQR(WANTT,WANTZ,N,ILO,IHI,H,LDH,W,ILOZ,IHIZ,
C    *                  Z,LDZ,INFO)
C
C  Purpose
C  =======
C
C  ZLAHQR is an auxiliary routine called by ZHSEQR to update the
C  eigenvalues and Schur decomposition already computed by ZHSEQR, by
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
C          It is assumed that H is already upper triangular in rows and
C          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
C          ZLAHQR works primarily with the Hessenberg submatrix in rows
C          and columns ILO to IHI, but applies transformations to all of
C          H if WANTT is .TRUE..
C          1 <= ILO <= max(1,IHI); IHI <= N.
C
C  H       (input/output) COMPLEX*16 array, dimension (LDH,N)
C          On entry, the upper Hessenberg matrix H.
C          On exit, if WANTT is .TRUE., H is upper triangular in rows
C          and columns ILO:IHI, with any 2-by-2 diagonal blocks in
C          standard form. If WANTT is .FALSE., the contents of H are
C          unspecified on exit.
C
C  LDH     (input) INTEGER
C          The leading dimension of the array H. LDH >= max(1,N).
C
C  W       (output) COMPLEX*16 array, dimension (N)
C          The computed eigenvalues ILO to IHI are stored in the
C          corresponding elements of W. If WANTT is .TRUE., the
C          eigenvalues are stored in the same order as on the diagonal
C          of the Schur form returned in H, with W(i) = H(i,i).
C
C  ILOZ    (input) INTEGER
C  IHIZ    (input) INTEGER
C          Specify the rows of Z to which transformations must be
C          applied if WANTZ is .TRUE..
C          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
C
C  Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)
C          If WANTZ is .TRUE., on entry Z must contain the current
C          matrix Z of transformations accumulated by ZHSEQR, and on
C          exit Z has been updated; transformations are applied only to
C          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
C          If WANTZ is .FALSE., Z is not referenced.
C
C  LDZ     (input) INTEGER
C          The leading dimension of the array Z. LDZ >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          > 0: ZLAHQR failed to compute all the eigenvalues ILO to IHI
C               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
C               elements i+1:ihi of W contain those eigenvalues which
C               have been successfully computed.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  RZERO, RONE, HALF
      PARAMETER         (RZERO=0.0D+0,RONE=1.0D+0,HALF=0.5D+0)
C     .. Scalar Arguments ..
      INTEGER           IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL           WANTT, WANTZ
C     .. Array Arguments ..
      COMPLEX*16        H(LDH,*), W(*), Z(LDZ,*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM, H11, H11S, H22, SUM, T, T1, TEMP, U, V2,
     *                  X, Y
      DOUBLE PRECISION  H10, H21, OVFL, RTEMP, S, SMLNUM, T2, TST1, ULP,
     *                  UNFL
      INTEGER           I, I1, I2, ITN, ITS, J, K, L, M, NH, NZ
      LOGICAL           DIVFLG
C     .. Local Arrays ..
      COMPLEX*16        V(2)
      DOUBLE PRECISION  RWORK(1)
C     .. External Functions ..
      COMPLEX*16        F06CLF
      DOUBLE PRECISION  F06BNF, F06UMF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          F06CLF, F06BNF, F06UMF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F08ASV, ZCOPY, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(CDUM) = ABS(DBLE(CDUM)) + ABS(DIMAG(CDUM))
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
      IF (ILO.EQ.IHI) THEN
         W(ILO) = H(ILO,ILO)
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
      OVFL = RONE/UNFL
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
C     IHI to ILO in steps of 1. Each iteration of the loop works
C     with the active submatrix in rows and columns L to I.
C     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
C     H(L,L-1) is negligible so that the matrix splits.
C
      I = IHI
   20 CONTINUE
      IF (I.LT.ILO) GO TO 260
C
C     Perform QR iterations on rows and columns ILO to I until a
C     submatrix of order 1 splits off at the bottom because a
C     subdiagonal element has become negligible.
C
      L = ILO
      DO 220 ITS = 0, ITN
C
C        Look for a single small subdiagonal element.
C
         DO 40 K = I, L + 1, -1
            TST1 = CABS1(H(K-1,K-1)) + CABS1(H(K,K))
            IF (TST1.EQ.RZERO) TST1 = F06UMF('1',I-L+1,H(L,L),LDH,RWORK)
            IF (ABS(DBLE(H(K,K-1))).LE.MAX(ULP*TST1,SMLNUM)) GO TO 60
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
C        Exit from loop if a submatrix of order 1 has split off.
C
         IF (L.GE.I) GO TO 240
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
            T = ABS(DBLE(H(I,I-1))) + ABS(DBLE(H(I-1,I-2)))
         ELSE
C
C           Wilkinson's shift.
C
            T = H(I,I)
            U = H(I-1,I)*DBLE(H(I,I-1))
            IF (U.NE.ZERO) THEN
               X = HALF*(H(I-1,I-1)-T)
               Y = SQRT(X*X+U)
               IF (DBLE(X)*DBLE(Y)+DIMAG(X)*DIMAG(Y).LT.RZERO) Y = -Y
               T = T - F06CLF(U,(X+Y),DIVFLG)
            END IF
         END IF
C
C        Look for two consecutive small subdiagonal elements.
C
         DO 80 M = I - 1, L, -1
C
C           Determine the effect of starting the single-shift QR
C           iteration at row M, and see if this would make H(M,M-1)
C           negligible.
C
            H11 = H(M,M)
            H22 = H(M+1,M+1)
            H11S = H11 - T
            H21 = H(M+1,M)
            S = CABS1(H11S) + ABS(H21)
            H11S = H11S/S
            H21 = H21/S
            V(1) = H11S
            V(2) = H21
            IF (M.EQ.L) GO TO 100
            H10 = H(M,M-1)
            TST1 = CABS1(H11S)*(CABS1(H11)+CABS1(H22))
            IF (ABS(H10*H21).LE.ULP*TST1) GO TO 100
   80    CONTINUE
  100    CONTINUE
C
C        Single-shift QR step
C
         DO 200 K = M, I - 1
C
C           The first iteration of this loop determines a reflection G
C           from the vector V and applies it from left and right to H,
C           thus creating a nonzero bulge below the subdiagonal.
C
C           Each subsequent iteration determines a reflection G to
C           restore the Hessenberg form in the (K-1)th column, and thus
C           chases the bulge one step toward the bottom of the active
C           submatrix.
C
C           V(2) is always real before the call to F08ASV, and hence
C           after the call T2 ( = T1*V(2) ) is also real.
C
            IF (K.GT.M) CALL ZCOPY(2,H(K,K-1),1,V,1)
            CALL F08ASV(2,V(1),V(2),1,T1)
            IF (K.GT.M) THEN
               H(K,K-1) = V(1)
               H(K+1,K-1) = ZERO
            END IF
            V2 = V(2)
            T2 = DBLE(T1*V2)
C
C           Apply G from the left to transform the rows of the matrix
C           in columns K to I2.
C
            DO 120 J = K, I2
               SUM = DCONJG(T1)*H(K,J) + T2*H(K+1,J)
               H(K,J) = H(K,J) - SUM
               H(K+1,J) = H(K+1,J) - SUM*V2
  120       CONTINUE
C
C           Apply G from the right to transform the columns of the
C           matrix in rows I1 to min(K+2,I).
C
            DO 140 J = I1, MIN(K+2,I)
               SUM = T1*H(J,K) + T2*H(J,K+1)
               H(J,K) = H(J,K) - SUM
               H(J,K+1) = H(J,K+1) - SUM*DCONJG(V2)
  140       CONTINUE
C
            IF (WANTZ) THEN
C
C              Accumulate transformations in the matrix Z
C
               DO 160 J = ILOZ, IHIZ
                  SUM = T1*Z(J,K) + T2*Z(J,K+1)
                  Z(J,K) = Z(J,K) - SUM
                  Z(J,K+1) = Z(J,K+1) - SUM*DCONJG(V2)
  160          CONTINUE
            END IF
C
            IF (K.EQ.M .AND. M.GT.L) THEN
C
C              If the QR step was started at row M > L because two
C              consecutive small subdiagonals were found, then extra
C              scaling must be performed to ensure that H(M,M-1) remains
C              real.
C
               TEMP = ONE - T1
               TEMP = TEMP/F06BNF(DBLE(TEMP),DIMAG(TEMP))
               H(M+1,M) = H(M+1,M)*DCONJG(TEMP)
               IF (M+2.LE.I) H(M+2,M+1) = H(M+2,M+1)*TEMP
               DO 180 J = M, I
                  IF (J.NE.M+1) THEN
                     IF (I2.GT.J) CALL ZSCAL(I2-J,TEMP,H(J,J+1),LDH)
                     CALL ZSCAL(J-I1,DCONJG(TEMP),H(I1,J),1)
                     IF (WANTZ) THEN
                        CALL ZSCAL(NZ,DCONJG(TEMP),Z(ILOZ,J),1)
                     END IF
                  END IF
  180          CONTINUE
            END IF
  200    CONTINUE
C
C        Ensure that H(I,I-1) is real.
C
         TEMP = H(I,I-1)
         IF (DIMAG(TEMP).NE.RZERO) THEN
            RTEMP = F06BNF(DBLE(TEMP),DIMAG(TEMP))
            H(I,I-1) = RTEMP
            TEMP = TEMP/RTEMP
            IF (I2.GT.I) CALL ZSCAL(I2-I,DCONJG(TEMP),H(I,I+1),LDH)
            CALL ZSCAL(I-I1,TEMP,H(I1,I),1)
            IF (WANTZ) THEN
               CALL ZSCAL(NZ,TEMP,Z(ILOZ,I),1)
            END IF
         END IF
C
  220 CONTINUE
C
C     Failure to converge in remaining number of iterations
C
      INFO = I
      RETURN
C
  240 CONTINUE
C
C     H(I,I-1) is negligible: one eigenvalue has converged.
C
      W(I) = H(I,I)
C
C     Decrement number of remaining iterations, and return to start of
C     the main loop with new value of I.
C
      ITN = ITN - ITS
      I = L - 1
      GO TO 20
C
  260 CONTINUE
      RETURN
C
C     End of F08PSZ (ZLAHQR)
C
      END
