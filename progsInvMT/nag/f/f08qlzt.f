      SUBROUTINE F08QLZ(LTRAN,LREAL,N,T,LDT,B,W,SCALE,X,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAQTR(LTRAN,LREAL,N,T,LDT,B,W,SCALE,X,WORK,
C    *                  INFO)
C
C  Purpose
C  =======
C
C  DLAQTR solves the real quasi-triangular system
C
C     op(T)*p = scale*c,               if LREAL = .TRUE.
C
C  or the complex quasi-triangular systems
C
C     op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.
C
C  in real arithmetic, where T is upper quasi-triangular.
C  If LREAL = .FALSE., then the first diagonal block of T
C  must be 1 by 1, B is the specially structured matrix
C
C                 B = [ b(1) b(2) ... b(n) ]
C                     [       w            ]
C                     [           w        ]
C                     [              .     ]
C                     [                 w  ]
C
C  op(A) = A or A', A' denotes the conjugate transpose of
C  matrix A.
C
C  On input, X = [ c ].  On output, X = [ p ].
C                [ d ]                  [ q ]
C
C  This subroutine is designed for the condition number estimation
C  in routine DTRSNA.
C
C  Arguments
C  =========
C
C  LTRAN   (input) LOGICAL
C          On entry, LTRAN specifies the option of conjugate
C          transpose:
C             = .FALSE.,    op(T+i*B) = T+i*B,
C             = .TRUE.,     op(T+i*B) = (T+i*B)'.
C
C  LREAL   (input) LOGICAL
C          On entry, LREAL specifies the input matrix
C          structure:
C             = .FALSE.,    the input is complex
C             = .TRUE.,     the input is real
C
C  N       (input) INTEGER
C          On entry, N specifies the order of T+i*B. N >= 0.
C
C  T       (input) DOUBLE PRECISION array, dimension(LDT,N)
C          On entry, T contains a matrix in Schur canonical form.
C          If LREAL = .FALSE., then the first diagonal block
C          of T must be 1 by 1.
C
C  LDT     (input) INTEGER
C          The leading dimension of the matrix T. LDT >= max(1,N).
C
C  B       (input) DOUBLE PRECISION array, dimension(N)
C          On entry, B contains the elements to form the matrix
C          B as described above.
C          If LREAL = .TRUE., B is not referenced.
C
C  W       (input) DOUBLE PRECISION
C          On entry, W is the diagonal element of the matrix B.
C          If LREAL = .TRUE., W is not referenced.
C
C  SCALE   (output) DOUBLE PRECISION
C          On exit, SCALE is the scale factor.
C
C  X       (input/output) DOUBLE PRECISION array, dimension(2*N)
C          On entry, X contains the right hand side of the system.
C          On exit, X is overwritten by the solution.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension(N)
C
C  INFO    (output) INTEGER
C          On exit, INFO is set to
C             0: successful exit.
C               1: the some diagonal 1 by 1 block has been perturbed by
C                  a small number SMIN to keep nonsingularity.
C               2: the some diagonal 2 by 2 block has been perturbed by
C                  a small number in F08QHX to keep nonsingularity.
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
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SCALE, W
      INTEGER           INFO, LDT, N
      LOGICAL           LREAL, LTRAN
C     .. Array Arguments ..
      DOUBLE PRECISION  B(*), T(LDT,*), WORK(*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGNUM, EPS, REC, SCALOC, SI, SMIN, SMINW,
     *                  SMLNUM, SR, TJJ, TMP, XJ, XMAX, XNORM, Z
      INTEGER           I, IERR, J, J1, J2, JNEXT, K, N1, N2
      LOGICAL           NOTRAN
C     .. Local Arrays ..
      DOUBLE PRECISION  D(2,2), V(2,2)
C     .. External Functions ..
      DOUBLE PRECISION  DASUM, DDOT, F06RAF, X02AJF, X02AMF
      INTEGER           IDAMAX, X02BHF
      EXTERNAL          DASUM, DDOT, F06RAF, X02AJF, X02AMF, IDAMAX,
     *                  X02BHF
C     .. External Subroutines ..
      EXTERNAL          A02ACF, DAXPY, DSCAL, F08QHX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
C     Do not test the input parameters for errors
C
      NOTRAN = .NOT. LTRAN
      INFO = 0
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Set constants to control overflow
C
      EPS = X02AJF()*X02BHF()
      SMLNUM = X02AMF()/EPS
      BIGNUM = ONE/SMLNUM
C
      XNORM = F06RAF('M',N,N,T,LDT,D)
      IF ( .NOT. LREAL) XNORM = MAX(XNORM,ABS(W),F06RAF('M',N,1,B,N,D))
      SMIN = MAX(SMLNUM,EPS*XNORM)
C
C     Compute 1-norm of each column of strictly upper triangular
C     part of T to control overflow in triangular solver.
C
      WORK(1) = ZERO
      DO 20 J = 2, N
         WORK(J) = DASUM(J-1,T(1,J),1)
   20 CONTINUE
C
      IF ( .NOT. LREAL) THEN
         DO 40 I = 2, N
            WORK(I) = WORK(I) + ABS(B(I))
   40    CONTINUE
      END IF
C
      N2 = 2*N
      N1 = N
      IF ( .NOT. LREAL) N1 = N2
      K = IDAMAX(N1,X,1)
      XMAX = ABS(X(K))
      SCALE = ONE
C
      IF (XMAX.GT.BIGNUM) THEN
         SCALE = BIGNUM/XMAX
         CALL DSCAL(N1,SCALE,X,1)
         XMAX = BIGNUM
      END IF
C
      IF (LREAL) THEN
C
         IF (NOTRAN) THEN
C
C           Solve T*p = scale*c
C
            JNEXT = N
            DO 60 J = N, 1, -1
               IF (J.GT.JNEXT) GO TO 60
               J1 = J
               J2 = J
               JNEXT = J - 1
               IF (J.GT.1 .AND. T(J,J-1).NE.ZERO) THEN
                  J1 = J - 1
                  J2 = J
                  JNEXT = J - 2
               END IF
C
               IF (J1.EQ.J2) THEN
C
C                 Meet 1 by 1 diagonal block
C
C                 Scale to avoid overflow when computing
C                     x(j) = b(j)/T(j,j)
C
                  XJ = ABS(X(J1))
                  TJJ = ABS(T(J1,J1))
                  TMP = T(J1,J1)
                  IF (TJJ.LT.SMIN) THEN
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  END IF
C
                  IF (XJ.EQ.ZERO) GO TO 60
C
                  IF (TJJ.LT.ONE) THEN
                     IF (XJ.GT.BIGNUM*TJJ) THEN
                        REC = ONE/XJ
                        CALL DSCAL(N,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X(J1) = X(J1)/TMP
                  XJ = ABS(X(J1))
C
C                 Scale x if necessary to avoid overflow when adding a
C                 multiple of column j1 of T.
C
                  IF (XJ.GT.ONE) THEN
                     REC = ONE/XJ
                     IF (WORK(J1).GT.(BIGNUM-XMAX)*REC) THEN
                        CALL DSCAL(N,REC,X,1)
                        SCALE = SCALE*REC
                     END IF
                  END IF
                  IF (J1.GT.1) THEN
                     CALL DAXPY(J1-1,-X(J1),T(1,J1),1,X,1)
                     K = IDAMAX(J1-1,X,1)
                     XMAX = ABS(X(K))
                  END IF
C
               ELSE
C
C                 Meet 2 by 2 diagonal block
C
C                 Call 2 by 2 linear system solve, to take
C                 care of possible overflow by scaling factor.
C
                  D(1,1) = X(J1)
                  D(2,1) = X(J2)
                  CALL F08QHX(.FALSE.,2,1,SMIN,ONE,T(J1,J1),LDT,ONE,ONE,
     *                        D,2,ZERO,ZERO,V,2,SCALOC,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 2
C
                  IF (SCALOC.NE.ONE) THEN
                     CALL DSCAL(N,SCALOC,X,1)
                     SCALE = SCALE*SCALOC
                  END IF
                  X(J1) = V(1,1)
                  X(J2) = V(2,1)
C
C                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
C                 to avoid overflow in updating right-hand side.
C
                  XJ = MAX(ABS(V(1,1)),ABS(V(2,1)))
                  IF (XJ.GT.ONE) THEN
                     REC = ONE/XJ
                     IF (MAX(WORK(J1),WORK(J2)).GT.(BIGNUM-XMAX)*REC)
     *                   THEN
                        CALL DSCAL(N,REC,X,1)
                        SCALE = SCALE*REC
                     END IF
                  END IF
C
C                 Update right-hand side
C
                  IF (J1.GT.1) THEN
                     CALL DAXPY(J1-1,-X(J1),T(1,J1),1,X,1)
                     CALL DAXPY(J1-1,-X(J2),T(1,J2),1,X,1)
                     K = IDAMAX(J1-1,X,1)
                     XMAX = ABS(X(K))
                  END IF
C
               END IF
C
   60       CONTINUE
C
         ELSE
C
C           Solve T'*p = scale*c
C
            JNEXT = 1
            DO 80 J = 1, N
               IF (J.LT.JNEXT) GO TO 80
               J1 = J
               J2 = J
               JNEXT = J + 1
               IF (J.LT.N .AND. T(J+1,J).NE.ZERO) THEN
                  J1 = J
                  J2 = J + 1
                  JNEXT = J + 2
               END IF
C
               IF (J1.EQ.J2) THEN
C
C                 1 by 1 diagonal block
C
C                 Scale if necessary to avoid overflow in forming the
C                 right-hand side entry by inner product.
C
                  XJ = ABS(X(J1))
                  IF (XMAX.GT.ONE) THEN
                     REC = ONE/XMAX
                     IF (WORK(J1).GT.(BIGNUM-XJ)*REC) THEN
                        CALL DSCAL(N,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
C
                  X(J1) = X(J1) - DDOT(J1-1,T(1,J1),1,X,1)
C
                  XJ = ABS(X(J1))
                  TJJ = ABS(T(J1,J1))
                  TMP = T(J1,J1)
                  IF (TJJ.LT.SMIN) THEN
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  END IF
C
                  IF (TJJ.LT.ONE) THEN
                     IF (XJ.GT.BIGNUM*TJJ) THEN
                        REC = ONE/XJ
                        CALL DSCAL(N,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X(J1) = X(J1)/TMP
                  XMAX = MAX(XMAX,ABS(X(J1)))
C
               ELSE
C
C                 2 by 2 diagonal block
C
C                 Scale if necessary to avoid overflow in forming the
C                 right-hand side entries by inner product.
C
                  XJ = MAX(ABS(X(J1)),ABS(X(J2)))
                  IF (XMAX.GT.ONE) THEN
                     REC = ONE/XMAX
                     IF (MAX(WORK(J2),WORK(J1)).GT.(BIGNUM-XJ)*REC) THEN
                        CALL DSCAL(N,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
C
                  D(1,1) = X(J1) - DDOT(J1-1,T(1,J1),1,X,1)
                  D(2,1) = X(J2) - DDOT(J1-1,T(1,J2),1,X,1)
C
                  CALL F08QHX(.TRUE.,2,1,SMIN,ONE,T(J1,J1),LDT,ONE,ONE,
     *                        D,2,ZERO,ZERO,V,2,SCALOC,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 2
C
                  IF (SCALOC.NE.ONE) THEN
                     CALL DSCAL(N,SCALOC,X,1)
                     SCALE = SCALE*SCALOC
                  END IF
                  X(J1) = V(1,1)
                  X(J2) = V(2,1)
                  XMAX = MAX(ABS(X(J1)),ABS(X(J2)),XMAX)
C
               END IF
   80       CONTINUE
         END IF
C
      ELSE
C
         SMINW = MAX(EPS*ABS(W),SMIN)
         IF (NOTRAN) THEN
C
C           Solve (T + iB)*(p+iq) = c+id
C
            JNEXT = N
            DO 140 J = N, 1, -1
               IF (J.GT.JNEXT) GO TO 140
               J1 = J
               J2 = J
               JNEXT = J - 1
               IF (J.GT.1 .AND. T(J,J-1).NE.ZERO) THEN
                  J1 = J - 1
                  J2 = J
                  JNEXT = J - 2
               END IF
C
               IF (J1.EQ.J2) THEN
C
C                 1 by 1 diagonal block
C
C                 Scale if necessary to avoid overflow in division
C
                  Z = W
                  IF (J1.EQ.1) Z = B(1)
                  XJ = ABS(X(J1)) + ABS(X(N+J1))
                  TJJ = ABS(T(J1,J1)) + ABS(Z)
                  TMP = T(J1,J1)
                  IF (TJJ.LT.SMINW) THEN
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  END IF
C
                  IF (XJ.EQ.ZERO) GO TO 140
C
                  IF (TJJ.LT.ONE) THEN
                     IF (XJ.GT.BIGNUM*TJJ) THEN
                        REC = ONE/XJ
                        CALL DSCAL(N2,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  CALL A02ACF(X(J1),X(N+J1),TMP,Z,SR,SI)
                  X(J1) = SR
                  X(N+J1) = SI
                  XJ = ABS(X(J1)) + ABS(X(N+J1))
C
C                 Scale x if necessary to avoid overflow when adding a
C                 multiple of column j1 of T.
C
                  IF (XJ.GT.ONE) THEN
                     REC = ONE/XJ
                     IF (WORK(J1).GT.(BIGNUM-XMAX)*REC) THEN
                        CALL DSCAL(N2,REC,X,1)
                        SCALE = SCALE*REC
                     END IF
                  END IF
C
                  IF (J1.GT.1) THEN
                     CALL DAXPY(J1-1,-X(J1),T(1,J1),1,X,1)
                     CALL DAXPY(J1-1,-X(N+J1),T(1,J1),1,X(N+1),1)
C
                     X(1) = X(1) + B(J1)*X(N+J1)
                     X(N+1) = X(N+1) - B(J1)*X(J1)
C
                     XMAX = ZERO
                     DO 100 K = 1, J1 - 1
                        XMAX = MAX(XMAX,ABS(X(K))+ABS(X(K+N)))
  100                CONTINUE
                  END IF
C
               ELSE
C
C                 Meet 2 by 2 diagonal block
C
                  D(1,1) = X(J1)
                  D(2,1) = X(J2)
                  D(1,2) = X(N+J1)
                  D(2,2) = X(N+J2)
                  CALL F08QHX(.FALSE.,2,2,SMINW,ONE,T(J1,J1),LDT,ONE,
     *                        ONE,D,2,ZERO,-W,V,2,SCALOC,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 2
C
                  IF (SCALOC.NE.ONE) THEN
                     CALL DSCAL(2*N,SCALOC,X,1)
                     SCALE = SCALOC*SCALE
                  END IF
                  X(J1) = V(1,1)
                  X(J2) = V(2,1)
                  X(N+J1) = V(1,2)
                  X(N+J2) = V(2,2)
C
C                 Scale X(J1), .... to avoid overflow in
C                 updating right hand side.
C
                  XJ = MAX(ABS(V(1,1))+ABS(V(1,2)),ABS(V(2,1))
     *                 +ABS(V(2,2)))
                  IF (XJ.GT.ONE) THEN
                     REC = ONE/XJ
                     IF (MAX(WORK(J1),WORK(J2)).GT.(BIGNUM-XMAX)*REC)
     *                   THEN
                        CALL DSCAL(N2,REC,X,1)
                        SCALE = SCALE*REC
                     END IF
                  END IF
C
C                 Update the right-hand side.
C
                  IF (J1.GT.1) THEN
                     CALL DAXPY(J1-1,-X(J1),T(1,J1),1,X,1)
                     CALL DAXPY(J1-1,-X(J2),T(1,J2),1,X,1)
C
                     CALL DAXPY(J1-1,-X(N+J1),T(1,J1),1,X(N+1),1)
                     CALL DAXPY(J1-1,-X(N+J2),T(1,J2),1,X(N+1),1)
C
                     X(1) = X(1) + B(J1)*X(N+J1) + B(J2)*X(N+J2)
                     X(N+1) = X(N+1) - B(J1)*X(J1) - B(J2)*X(J2)
C
                     XMAX = ZERO
                     DO 120 K = 1, J1 - 1
                        XMAX = MAX(ABS(X(K))+ABS(X(K+N)),XMAX)
  120                CONTINUE
                  END IF
C
               END IF
  140       CONTINUE
C
         ELSE
C
C           Solve (T + iB)'*(p+iq) = c+id
C
            JNEXT = 1
            DO 160 J = 1, N
               IF (J.LT.JNEXT) GO TO 160
               J1 = J
               J2 = J
               JNEXT = J + 1
               IF (J.LT.N .AND. T(J+1,J).NE.ZERO) THEN
                  J1 = J
                  J2 = J + 1
                  JNEXT = J + 2
               END IF
C
               IF (J1.EQ.J2) THEN
C
C                 1 by 1 diagonal block
C
C                 Scale if necessary to avoid overflow in forming the
C                 right-hand side entry by inner product.
C
                  XJ = ABS(X(J1)) + ABS(X(J1+N))
                  IF (XMAX.GT.ONE) THEN
                     REC = ONE/XMAX
                     IF (WORK(J1).GT.(BIGNUM-XJ)*REC) THEN
                        CALL DSCAL(N2,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
C
                  X(J1) = X(J1) - DDOT(J1-1,T(1,J1),1,X,1)
                  X(N+J1) = X(N+J1) - DDOT(J1-1,T(1,J1),1,X(N+1),1)
                  IF (J1.GT.1) THEN
                     X(J1) = X(J1) - B(J1)*X(N+1)
                     X(N+J1) = X(N+J1) + B(J1)*X(1)
                  END IF
                  XJ = ABS(X(J1)) + ABS(X(J1+N))
C
                  Z = W
                  IF (J1.EQ.1) Z = B(1)
C
C                 Scale if necessary to avoid overflow in
C                 complex division
C
                  TJJ = ABS(T(J1,J1)) + ABS(Z)
                  TMP = T(J1,J1)
                  IF (TJJ.LT.SMINW) THEN
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  END IF
C
                  IF (TJJ.LT.ONE) THEN
                     IF (XJ.GT.BIGNUM*TJJ) THEN
                        REC = ONE/XJ
                        CALL DSCAL(N2,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  CALL A02ACF(X(J1),X(N+J1),TMP,-Z,SR,SI)
                  X(J1) = SR
                  X(J1+N) = SI
                  XMAX = MAX(ABS(X(J1))+ABS(X(J1+N)),XMAX)
C
               ELSE
C
C                 2 by 2 diagonal block
C
C                 Scale if necessary to avoid overflow in forming the
C                 right-hand side entry by inner product.
C
                  XJ = MAX(ABS(X(J1))+ABS(X(N+J1)),ABS(X(J2))
     *                 +ABS(X(N+J2)))
                  IF (XMAX.GT.ONE) THEN
                     REC = ONE/XMAX
                     IF (MAX(WORK(J1),WORK(J2)).GT.(BIGNUM-XJ)/XMAX)
     *                   THEN
                        CALL DSCAL(N2,REC,X,1)
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
C
                  D(1,1) = X(J1) - DDOT(J1-1,T(1,J1),1,X,1)
                  D(2,1) = X(J2) - DDOT(J1-1,T(1,J2),1,X,1)
                  D(1,2) = X(N+J1) - DDOT(J1-1,T(1,J1),1,X(N+1),1)
                  D(2,2) = X(N+J2) - DDOT(J1-1,T(1,J2),1,X(N+1),1)
                  D(1,1) = D(1,1) - B(J1)*X(N+1)
                  D(2,1) = D(2,1) - B(J2)*X(N+1)
                  D(1,2) = D(1,2) + B(J1)*X(1)
                  D(2,2) = D(2,2) + B(J2)*X(1)
C
                  CALL F08QHX(.TRUE.,2,2,SMINW,ONE,T(J1,J1),LDT,ONE,ONE,
     *                        D,2,ZERO,W,V,2,SCALOC,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 2
C
                  IF (SCALOC.NE.ONE) THEN
                     CALL DSCAL(N2,SCALOC,X,1)
                     SCALE = SCALOC*SCALE
                  END IF
                  X(J1) = V(1,1)
                  X(J2) = V(2,1)
                  X(N+J1) = V(1,2)
                  X(N+J2) = V(2,2)
                  XMAX = MAX(ABS(X(J1))+ABS(X(N+J1)),ABS(X(J2))
     *                   +ABS(X(N+J2)),XMAX)
C
               END IF
C
  160       CONTINUE
C
         END IF
C
      END IF
C
      RETURN
C
C     End of F08QLZ (DLAQTR)
C
      END
