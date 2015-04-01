      SUBROUTINE F08QTF(COMPQ,N,T,LDT,Q,LDQ,IFST,ILST,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZTREXC(COMPQ,N,T,LDT,Q,LDQ,IFST,ILST,INFO)
C
C  Purpose
C  =======
C
C  ZTREXC reorders the Schur factorization of a complex matrix
C  A = Q*T*Q', so that the diagonal element of T with row index IFST is
C  moved to row ILST.
C
C  The Schur form T is reordered by a unitary similarity transformation
C  Z'*T*Z, and optionally the matrix Q of Schur vectors is updated by
C  postmultplying it with Z.
C
C  Arguments
C  =========
C
C  COMPQ   (input) CHARACTER*1
C          = 'V': update the matrix Q of Schur vectors;
C          = 'N': do not update Q.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input/output) COMPLEX*16 array, dimension (LDT,N)
C          On entry, the upper triangular matrix T.
C          On exit, T is overwritten by the reordered matrix T.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)
C          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
C          On exit, if COMPQ = 'V', Q has been postmultiplied by the
C          unitary transformation matrix Z which reorders T.
C          If COMPQ = 'N', Q is not referenced.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q.
C          LDQ >= 1; and if COMPQ = 'V', LDQ >= max(1,N).
C
C  IFST    (input) INTEGER
C  ILST    (input) INTEGER
C          Specify the re-ordering of the diagonal elements of T:
C          The element with row index IFST is moved to row ILST by a
C          sequence of transpositions between adjacent elements.
C          1 <= IFST <= N; 1 <= ILST <= N.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  RZERO
      PARAMETER         (RZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFST, ILST, INFO, LDQ, LDT, N
      CHARACTER         COMPQ
C     .. Array Arguments ..
      COMPLEX*16        Q(LDQ,*), T(LDT,*)
C     .. Local Scalars ..
      COMPLEX*16        SN, T11, T22, TEMP
      DOUBLE PRECISION  CS, QMAX, RTEMP
      INTEGER           I, II, J, K, M1, M2, M3
      LOGICAL           WANTQ
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07MRW, F08HSW, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCONJG, MAX
C     .. Executable Statements ..
C
C     Decode and test the input parameters.
C
      INFO = 0
      WANTQ = (COMPQ.EQ.'V' .OR. COMPQ.EQ.'v')
      IF ( .NOT. (COMPQ.EQ.'N' .OR. COMPQ.EQ.'n') .AND. .NOT. WANTQ)
     *    THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDT.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LDQ.LT.1 .OR. (WANTQ .AND. LDQ.LT.MAX(1,N))) THEN
         INFO = -6
      ELSE IF (IFST.LT.1 .OR. IFST.GT.N) THEN
         INFO = -7
      ELSE IF (ILST.LT.1 .OR. ILST.GT.N) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QTF/ZTREXC',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.1 .OR. IFST.EQ.ILST) RETURN
C
      IF (IFST.LT.ILST) THEN
C
C        Move the IFST-th diagonal element forward down the diagonal.
C
         M1 = 0
         M2 = -1
         M3 = 1
      ELSE
C
C        Move the IFST-th diagonal element backward up the diagonal.
C
         M1 = -1
         M2 = 0
         M3 = -1
      END IF
C
      DO 20 K = IFST + M1, ILST + M2, M3
C
C        Interchange the k-th and (k+1)-th diagonal elements.
C
         T11 = T(K,K)
         T22 = T(K+1,K+1)
C
C        Determine the transformation to perform the interchange.
C
         CALL F08HSW(T(K,K+1),T22-T11,CS,SN,TEMP)
C
C        Apply transformation to the matrix T.
C
         IF (K+2.LE.N) CALL F07MRW(N-K-1,T(K,K+2),LDT,T(K+1,K+2),LDT,CS,
     *                             SN)
         CALL F07MRW(K-1,T(1,K),1,T(1,K+1),1,CS,DCONJG(SN))
C
         T(K,K) = T22
         T(K+1,K+1) = T11
C
         IF (WANTQ) THEN
C
C           Accumulate transformation in the matrix Q.
C
            CALL F07MRW(N,Q(1,K),1,Q(1,K+1),1,CS,DCONJG(SN))
         END IF
C
   20 CONTINUE
C
      IF (WANTQ) THEN
C
C        Normalize Schur vectors so that element of largest absolute
C        value is real and positive.
C
         DO 60 J = 1, N
            QMAX = RZERO
            DO 40 I = 1, N
               RTEMP = ABS(Q(I,J))
               IF (RTEMP.GT.QMAX) THEN
                  II = I
                  QMAX = RTEMP
               END IF
   40       CONTINUE
            TEMP = Q(II,J)/QMAX
            CALL ZSCAL(N,DCONJG(TEMP),Q(1,J),1)
            Q(II,J) = QMAX
            CALL ZSCAL(J-1,DCONJG(TEMP),T(1,J),1)
            IF (J.LT.N) CALL ZSCAL(N-J,TEMP,T(J,J+1),LDT)
   60    CONTINUE
      END IF
      RETURN
C
C     End of F08QTF (ZTREXC)
C
      END
