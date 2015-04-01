      SUBROUTINE F08ASY(SIDE,TRANS,DIRECT,STOREV,M,N,K,V,LDV,T,LDT,C,
     *                  LDC,WORK,LDWORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARFB(SIDE,TRANS,DIRECT,STOREV,M,N,K,V,LDV,T,
C    *                  LDT,C,LDC,WORK,LDWORK)
C
C  Purpose
C  =======
C
C  ZLARFB applies a complex block reflector H or its transpose H' to a
C  complex m by n matrix C, from either the left or the right.
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': apply H or H' from the Left
C          = 'R': apply H or H' from the Right
C
C  TRANS   (input) CHARACTER*1
C          = 'N': apply H (No transpose)
C          = 'C': apply H' (Conjugate transpose)
C
C  DIRECT  (input) CHARACTER*1
C          Indicates how H is formed from a product of elementary
C          reflectors
C          = 'F': H = H(1) H(2) . . . H(k) (Forward)
C          = 'B': H = H(k) . . . H(2) H(1) (Backward)
C
C  STOREV  (input) CHARACTER*1
C          Indicates how the vectors which define the elementary
C          reflectors are stored:
C          = 'C': Columnwise
C          = 'R': Rowwise
C
C  M       (input) INTEGER
C          The number of rows of the matrix C.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C.
C
C  K       (input) INTEGER
C          The order of the matrix T (= the number of elementary
C          reflectors whose product defines the block reflector).
C
C  V       (input) COMPLEX*16 array, dimension
C                                (LDV,K) if STOREV = 'C'
C                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
C                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
C          The matrix V. See further details.
C
C  LDV     (input) INTEGER
C          The leading dimension of the array V.
C          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
C          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
C          if STOREV = 'R', LDV >= K.
C
C  T       (input) COMPLEX*16 array, dimension (LDT,K)
C          The triangular k-by-k matrix T in the representation of the
C          block reflector.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= K.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the m-by-n matrix C.
C          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDA >= max(1,M).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)
C
C  LDWORK  (input) INTEGER
C          The leading dimension of the array WORK.
C          If SIDE = 'L', LDWORK >= max(1,N);
C          if SIDE = 'R', LDWORK >= max(1,M).
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K, LDC, LDT, LDV, LDWORK, M, N
      CHARACTER         DIRECT, SIDE, STOREV, TRANS
C     .. Array Arguments ..
      COMPLEX*16        C(LDC,*), T(LDT,*), V(LDV,*), WORK(LDWORK,*)
C     .. Local Scalars ..
      INTEGER           I, J
      CHARACTER         TRANST
C     .. External Subroutines ..
      EXTERNAL          F07FRY, ZCOPY, ZGEMM, ZTRMM
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF (M.LE.0 .OR. N.LE.0) RETURN
C
      IF ((TRANS.EQ.'N' .OR. TRANS.EQ.'n')) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
C
      IF ((STOREV.EQ.'C' .OR. STOREV.EQ.'c')) THEN
C
         IF ((DIRECT.EQ.'F' .OR. DIRECT.EQ.'f')) THEN
C
C           Let  V =  ( V1 )    (first K rows)
C                     ( V2 )
C           where  V1  is unit lower triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
C
C              W := C1'
C
               DO 20 J = 1, K
                  CALL ZCOPY(N,C(J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
   20          CONTINUE
C
C              W := W * V1
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',N,K,ONE,
     *                    V,LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C2'*V2
C
                  CALL ZGEMM('Conjugate transpose','No transpose',N,K,
     *                       M-K,ONE,C(K+1,1),LDC,V(K+1,1),LDV,ONE,WORK,
     *                       LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Upper',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V * W'
C
               IF (M.GT.K) THEN
C
C                 C2 := C2 - V2 * W'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M-K,N,
     *                       K,-ONE,V(K+1,1),LDV,WORK,LDWORK,ONE,
     *                       C(K+1,1),LDC)
               END IF
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    N,K,ONE,V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W'
C
               DO 60 J = 1, K
                  DO 40 I = 1, N
                     C(J,I) = C(J,I) - DCONJG(WORK(I,J))
   40             CONTINUE
   60          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C1
C
               DO 80 J = 1, K
                  CALL ZCOPY(M,C(1,J),1,WORK(1,J),1)
   80          CONTINUE
C
C              W := W * V1
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',M,K,ONE,
     *                    V,LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C2 * V2
C
                  CALL ZGEMM('No transpose','No transpose',M,K,N-K,ONE,
     *                       C(1,K+1),LDC,V(K+1,1),LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Upper',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V'
C
               IF (N.GT.K) THEN
C
C                 C2 := C2 - W * V2'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,N-K,
     *                       K,-ONE,WORK,LDWORK,V(K+1,1),LDV,ONE,
     *                       C(1,K+1),LDC)
               END IF
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    M,K,ONE,V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W
C
               DO 120 J = 1, K
                  DO 100 I = 1, M
                     C(I,J) = C(I,J) - WORK(I,J)
  100             CONTINUE
  120          CONTINUE
            END IF
C
         ELSE
C
C           Let  V =  ( V1 )
C                     ( V2 )    (last K rows)
C           where  V2  is unit upper triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
C
C              W := C2'
C
               DO 140 J = 1, K
                  CALL ZCOPY(N,C(M-K+J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
  140          CONTINUE
C
C              W := W * V2
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',N,K,ONE,
     *                    V(M-K+1,1),LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C1'*V1
C
                  CALL ZGEMM('Conjugate transpose','No transpose',N,K,
     *                       M-K,ONE,C,LDC,V,LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Lower',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V * W'
C
               IF (M.GT.K) THEN
C
C                 C1 := C1 - V1 * W'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M-K,N,
     *                       K,-ONE,V,LDV,WORK,LDWORK,ONE,C,LDC)
               END IF
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    N,K,ONE,V(M-K+1,1),LDV,WORK,LDWORK)
C
C              C2 := C2 - W'
C
               DO 180 J = 1, K
                  DO 160 I = 1, N
                     C(M-K+J,I) = C(M-K+J,I) - DCONJG(WORK(I,J))
  160             CONTINUE
  180          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C2
C
               DO 200 J = 1, K
                  CALL ZCOPY(M,C(1,N-K+J),1,WORK(1,J),1)
  200          CONTINUE
C
C              W := W * V2
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',M,K,ONE,
     *                    V(N-K+1,1),LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C1 * V1
C
                  CALL ZGEMM('No transpose','No transpose',M,K,N-K,ONE,
     *                       C,LDC,V,LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Lower',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V'
C
               IF (N.GT.K) THEN
C
C                 C1 := C1 - W * V1'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,N-K,
     *                       K,-ONE,WORK,LDWORK,V,LDV,ONE,C,LDC)
               END IF
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    M,K,ONE,V(N-K+1,1),LDV,WORK,LDWORK)
C
C              C2 := C2 - W
C
               DO 240 J = 1, K
                  DO 220 I = 1, M
                     C(I,N-K+J) = C(I,N-K+J) - WORK(I,J)
  220             CONTINUE
  240          CONTINUE
            END IF
         END IF
C
      ELSE IF ((STOREV.EQ.'R' .OR. STOREV.EQ.'r')) THEN
C
         IF ((DIRECT.EQ.'F' .OR. DIRECT.EQ.'f')) THEN
C
C           Let  V =  ( V1  V2 )    (V1: first K columns)
C           where  V1  is unit upper triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
C
C              W := C1'
C
               DO 260 J = 1, K
                  CALL ZCOPY(N,C(J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
  260          CONTINUE
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    N,K,ONE,V,LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C2'*V2'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',N,K,M-K,ONE,C(K+1,1),
     *                       LDC,V(1,K+1),LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Upper',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V' * W'
C
               IF (M.GT.K) THEN
C
C                 C2 := C2 - V2' * W'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',M-K,N,K,-ONE,V(1,K+1)
     *                       ,LDV,WORK,LDWORK,ONE,C(K+1,1),LDC)
               END IF
C
C              W := W * V1
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',N,K,ONE,
     *                    V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W'
C
               DO 300 J = 1, K
                  DO 280 I = 1, N
                     C(J,I) = C(J,I) - DCONJG(WORK(I,J))
  280             CONTINUE
  300          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
C
C              W := C1
C
               DO 320 J = 1, K
                  CALL ZCOPY(M,C(1,J),1,WORK(1,J),1)
  320          CONTINUE
C
C              W := W * V1'
C
               CALL ZTRMM('Right','Upper','Conjugate transpose','Unit',
     *                    M,K,ONE,V,LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C2 * V2'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,K,
     *                       N-K,ONE,C(1,K+1),LDC,V(1,K+1),LDV,ONE,WORK,
     *                       LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Upper',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V
C
               IF (N.GT.K) THEN
C
C                 C2 := C2 - W * V2
C
                  CALL ZGEMM('No transpose','No transpose',M,N-K,K,-ONE,
     *                       WORK,LDWORK,V(1,K+1),LDV,ONE,C(1,K+1),LDC)
               END IF
C
C              W := W * V1
C
               CALL ZTRMM('Right','Upper','No transpose','Unit',M,K,ONE,
     *                    V,LDV,WORK,LDWORK)
C
C              C1 := C1 - W
C
               DO 360 J = 1, K
                  DO 340 I = 1, M
                     C(I,J) = C(I,J) - WORK(I,J)
  340             CONTINUE
  360          CONTINUE
C
            END IF
C
         ELSE
C
C           Let  V =  ( V1  V2 )    (V2: last K columns)
C           where  V2  is unit lower triangular.
C
            IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C              Form  H * C  or  H' * C  where  C = ( C1 )
C                                                  ( C2 )
C
C              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
C
C              W := C2'
C
               DO 380 J = 1, K
                  CALL ZCOPY(N,C(M-K+J,1),LDC,WORK(1,J),1)
                  CALL F07FRY(N,WORK(1,J),1)
  380          CONTINUE
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    N,K,ONE,V(1,M-K+1),LDV,WORK,LDWORK)
               IF (M.GT.K) THEN
C
C                 W := W + C1'*V1'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',N,K,M-K,ONE,C,LDC,V,
     *                       LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T'  or  W * T
C
               CALL ZTRMM('Right','Lower',TRANST,'Non-unit',N,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - V' * W'
C
               IF (M.GT.K) THEN
C
C                 C1 := C1 - V1' * W'
C
                  CALL ZGEMM('Conjugate transpose',
     *                       'Conjugate transpose',M-K,N,K,-ONE,V,LDV,
     *                       WORK,LDWORK,ONE,C,LDC)
               END IF
C
C              W := W * V2
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',N,K,ONE,
     *                    V(1,M-K+1),LDV,WORK,LDWORK)
C
C              C2 := C2 - W'
C
               DO 420 J = 1, K
                  DO 400 I = 1, N
                     C(M-K+J,I) = C(M-K+J,I) - DCONJG(WORK(I,J))
  400             CONTINUE
  420          CONTINUE
C
            ELSE IF ((SIDE.EQ.'R' .OR. SIDE.EQ.'r')) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
C
C              W := C2
C
               DO 440 J = 1, K
                  CALL ZCOPY(M,C(1,N-K+J),1,WORK(1,J),1)
  440          CONTINUE
C
C              W := W * V2'
C
               CALL ZTRMM('Right','Lower','Conjugate transpose','Unit',
     *                    M,K,ONE,V(1,N-K+1),LDV,WORK,LDWORK)
               IF (N.GT.K) THEN
C
C                 W := W + C1 * V1'
C
                  CALL ZGEMM('No transpose','Conjugate transpose',M,K,
     *                       N-K,ONE,C,LDC,V,LDV,ONE,WORK,LDWORK)
               END IF
C
C              W := W * T  or  W * T'
C
               CALL ZTRMM('Right','Lower',TRANS,'Non-unit',M,K,ONE,T,
     *                    LDT,WORK,LDWORK)
C
C              C := C - W * V
C
               IF (N.GT.K) THEN
C
C                 C1 := C1 - W * V1
C
                  CALL ZGEMM('No transpose','No transpose',M,N-K,K,-ONE,
     *                       WORK,LDWORK,V,LDV,ONE,C,LDC)
               END IF
C
C              W := W * V2
C
               CALL ZTRMM('Right','Lower','No transpose','Unit',M,K,ONE,
     *                    V(1,N-K+1),LDV,WORK,LDWORK)
C
C              C1 := C1 - W
C
               DO 480 J = 1, K
                  DO 460 I = 1, M
                     C(I,N-K+J) = C(I,N-K+J) - WORK(I,J)
  460             CONTINUE
  480          CONTINUE
C
            END IF
C
         END IF
      END IF
C
      RETURN
C
C     End of F08ASY (ZLARFB)
C
      END
