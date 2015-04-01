      SUBROUTINE F08QFZ(WANTQ,N,T,LDT,Q,LDQ,J1,N1,N2,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAEXC(WANTQ,N,T,LDT,Q,LDQ,J1,N1,N2,WORK,INFO)
C
C  Purpose
C  =======
C
C  DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
C  an upper quasi-triangular matrix T by an orthogonal similarity
C  transformation.
C
C  T must be in Schur canonical form, that is, block upper triangular
C  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
C  has its diagonal elemnts equal and its off-diagonal elements of
C  opposite sign.
C
C  Arguments
C  =========
C
C  WANTQ   (input) LOGICAL
C          = .TRUE. : accumulate the transformation in the matrix Q;
C          = .FALSE.: do not accumulate the transformation.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C          On entry, the upper quasi-triangular matrix T, in Schur
C          canonical form.
C          On exit, the updated matrix T, again in Schur canonical form.
C
C  LDT     (input)  INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
C          On exit, if WANTQ is .TRUE., the updated matrix Q.
C          If WANTQ is .FALSE., Q is not referenced.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q.
C          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
C
C  J1      (input) INTEGER
C          The index of the first row of the first block T11.
C
C  N1      (input) INTEGER
C          The order of the first block T11. N1 = 0, 1 or 2.
C
C  N2      (input) INTEGER
C          The order of the second block T22. N2 = 0, 1 or 2.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          = 1: the transformed matrix T would be too far from Schur
C               form; the blocks are not swapped and T and Q are
C               unchanged.
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
      DOUBLE PRECISION  TEN
      PARAMETER         (TEN=1.0D+1)
      INTEGER           LDD, LDX
      PARAMETER         (LDD=4,LDX=2)
C     .. Scalar Arguments ..
      INTEGER           INFO, J1, LDQ, LDT, N, N1, N2
      LOGICAL           WANTQ
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,*), T(LDT,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22,
     *                  T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2,
     *                  WR1, WR2, XNORM
      INTEGER           IERR, J2, J3, J4, K, ND
C     .. Local Arrays ..
      DOUBLE PRECISION  D(LDD,4), U(3), U1(3), U2(3), X(LDX,2)
C     .. External Functions ..
      DOUBLE PRECISION  F06RAF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          F06RAF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          DROT, F06QFF, F08AEV, F08HEW, F08PEW, F08PEY,
     *                  F08QHY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0) RETURN
      IF (J1+N1.GT.N) RETURN
C
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
C
      IF (N1.EQ.1 .AND. N2.EQ.1) THEN
C
C        Swap two 1-by-1 blocks.
C
         T11 = T(J1,J1)
         T22 = T(J2,J2)
C
C        Determine the transformation to perform the interchange.
C
         CALL F08HEW(T(J1,J2),T22-T11,CS,SN,TEMP)
C
C        Apply transformation to the matrix T.
C
         IF (J3.LE.N) CALL DROT(N-J1-1,T(J1,J3),LDT,T(J2,J3),LDT,CS,SN)
         CALL DROT(J1-1,T(1,J1),1,T(1,J2),1,CS,SN)
C
         T(J1,J1) = T22
         T(J2,J2) = T11
C
         IF (WANTQ) THEN
C
C           Accumulate transformation in the matrix Q.
C
            CALL DROT(N,Q(1,J1),1,Q(1,J2),1,CS,SN)
         END IF
C
      ELSE
C
C        Swapping involves at least one 2-by-2 block.
C
C        Copy the diagonal block of order N1+N2 to the local array D
C        and compute its norm.
C
         ND = N1 + N2
         CALL F06QFF('General',ND,ND,T(J1,J1),LDT,D,LDD)
         DNORM = F06RAF('Max',ND,ND,D,LDD,WORK)
C
C        Compute machine-dependent threshold for test for accepting
C        swap.
C
         EPS = X02AJF()*X02BHF()
         SMLNUM = X02AMF()/EPS
         THRESH = MAX(TEN*EPS*DNORM,SMLNUM)
C
C        Solve T11*X - X*T22 = scale*T12 for X.
C
         CALL F08QHY(.FALSE.,.FALSE.,-1,N1,N2,D,LDD,D(N1+1,N1+1),LDD,
     *               D(1,N1+1),LDD,SCALE,X,LDX,XNORM,IERR)
C
C        Swap the adjacent diagonal blocks.
C
         K = N1 + N1 + N2 - 3
         GO TO (20,40,60) K
C
   20    CONTINUE
C
C        N1 = 1, N2 = 2: generate elementary reflector H so that:
C
C        ( scale, X11, X12 ) H = ( 0, 0, * )
C
         U(1) = SCALE
         U(2) = X(1,1)
         U(3) = X(1,2)
         CALL F08AEV(3,U(3),U,1,TAU)
         U(3) = ONE
         T11 = T(J1,J1)
C
C        Perform swap provisionally on diagonal block in D.
C
         CALL F08PEW('L',3,3,U,TAU,D,LDD,WORK)
         CALL F08PEW('R',3,3,U,TAU,D,LDD,WORK)
C
C        Test whether to reject swap.
C
         IF (MAX(ABS(D(3,1)),ABS(D(3,2)),ABS(D(3,3)-T11)).GT.THRESH)
     *       GO TO 100
C
C        Accept swap: apply transformation to the entire matrix T.
C
         CALL F08PEW('L',3,N-J1+1,U,TAU,T(J1,J1),LDT,WORK)
         CALL F08PEW('R',J2,3,U,TAU,T(1,J1),LDT,WORK)
C
         T(J3,J1) = ZERO
         T(J3,J2) = ZERO
         T(J3,J3) = T11
C
         IF (WANTQ) THEN
C
C           Accumulate transformation in the matrix Q.
C
            CALL F08PEW('R',N,3,U,TAU,Q(1,J1),LDQ,WORK)
         END IF
         GO TO 80
C
   40    CONTINUE
C
C        N1 = 2, N2 = 1: generate elementary reflector H so that:
C
C        H (  -X11 ) = ( * )
C          (  -X21 ) = ( 0 )
C          ( scale ) = ( 0 )
C
         U(1) = -X(1,1)
         U(2) = -X(2,1)
         U(3) = SCALE
         CALL F08AEV(3,U(1),U(2),1,TAU)
         U(1) = ONE
         T33 = T(J3,J3)
C
C        Perform swap provisionally on diagonal block in D.
C
         CALL F08PEW('L',3,3,U,TAU,D,LDD,WORK)
         CALL F08PEW('R',3,3,U,TAU,D,LDD,WORK)
C
C        Test whether to reject swap.
C
         IF (MAX(ABS(D(2,1)),ABS(D(3,1)),ABS(D(1,1)-T33)).GT.THRESH)
     *       GO TO 100
C
C        Accept swap: apply transformation to the entire matrix T.
C
         CALL F08PEW('R',J3,3,U,TAU,T(1,J1),LDT,WORK)
         CALL F08PEW('L',3,N-J1,U,TAU,T(J1,J2),LDT,WORK)
C
         T(J1,J1) = T33
         T(J2,J1) = ZERO
         T(J3,J1) = ZERO
C
         IF (WANTQ) THEN
C
C           Accumulate transformation in the matrix Q.
C
            CALL F08PEW('R',N,3,U,TAU,Q(1,J1),LDQ,WORK)
         END IF
         GO TO 80
C
   60    CONTINUE
C
C        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
C        that:
C
C        H(2) H(1) (  -X11  -X12 ) = (  *  * )
C                  (  -X21  -X22 )   (  0  * )
C                  ( scale    0  )   (  0  0 )
C                  (    0  scale )   (  0  0 )
C
         U1(1) = -X(1,1)
         U1(2) = -X(2,1)
         U1(3) = SCALE
         CALL F08AEV(3,U1(1),U1(2),1,TAU1)
         U1(1) = ONE
C
         TEMP = -TAU1*(X(1,2)+U1(2)*X(2,2))
         U2(1) = -TEMP*U1(2) - X(2,2)
         U2(2) = -TEMP*U1(3)
         U2(3) = SCALE
         CALL F08AEV(3,U2(1),U2(2),1,TAU2)
         U2(1) = ONE
C
C        Perform swap provisionally on diagonal block in D.
C
         CALL F08PEW('L',3,4,U1,TAU1,D,LDD,WORK)
         CALL F08PEW('R',4,3,U1,TAU1,D,LDD,WORK)
         CALL F08PEW('L',3,4,U2,TAU2,D(2,1),LDD,WORK)
         CALL F08PEW('R',4,3,U2,TAU2,D(1,2),LDD,WORK)
C
C        Test whether to reject swap.
C
         IF (MAX(ABS(D(3,1)),ABS(D(3,2)),ABS(D(4,1)),ABS(D(4,2)))
     *       .GT.THRESH) GO TO 100
C
C        Accept swap: apply transformation to the entire matrix T.
C
         CALL F08PEW('L',3,N-J1+1,U1,TAU1,T(J1,J1),LDT,WORK)
         CALL F08PEW('R',J4,3,U1,TAU1,T(1,J1),LDT,WORK)
         CALL F08PEW('L',3,N-J1+1,U2,TAU2,T(J2,J1),LDT,WORK)
         CALL F08PEW('R',J4,3,U2,TAU2,T(1,J2),LDT,WORK)
C
         T(J3,J1) = ZERO
         T(J3,J2) = ZERO
         T(J4,J1) = ZERO
         T(J4,J2) = ZERO
C
         IF (WANTQ) THEN
C
C           Accumulate transformation in the matrix Q.
C
            CALL F08PEW('R',N,3,U1,TAU1,Q(1,J1),LDQ,WORK)
            CALL F08PEW('R',N,3,U2,TAU2,Q(1,J2),LDQ,WORK)
         END IF
C
   80    CONTINUE
C
         IF (N2.EQ.2) THEN
C
C           Standardize new 2-by-2 block T11
C
            CALL F08PEY(T(J1,J1),T(J1,J2),T(J2,J1),T(J2,J2),WR1,WI1,WR2,
     *                  WI2,CS,SN)
            CALL DROT(N-J1-1,T(J1,J1+2),LDT,T(J2,J1+2),LDT,CS,SN)
            CALL DROT(J1-1,T(1,J1),1,T(1,J2),1,CS,SN)
            IF (WANTQ) CALL DROT(N,Q(1,J1),1,Q(1,J2),1,CS,SN)
         END IF
C
         IF (N1.EQ.2) THEN
C
C           Standardize new 2-by-2 block T22
C
            J3 = J1 + N2
            J4 = J3 + 1
            CALL F08PEY(T(J3,J3),T(J3,J4),T(J4,J3),T(J4,J4),WR1,WI1,WR2,
     *                  WI2,CS,SN)
            IF (J3+2.LE.N) CALL DROT(N-J3-1,T(J3,J3+2),LDT,T(J4,J3+2),
     *                               LDT,CS,SN)
            CALL DROT(J3-1,T(1,J3),1,T(1,J4),1,CS,SN)
            IF (WANTQ) CALL DROT(N,Q(1,J3),1,Q(1,J4),1,CS,SN)
         END IF
C
      END IF
      RETURN
C
C     Exit with INFO = 1 if swap was rejected.
C
  100 CONTINUE
      INFO = 1
      RETURN
C
C     End of F08QFZ (DLAEXC)
C
      END
