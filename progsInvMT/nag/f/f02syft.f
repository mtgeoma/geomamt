      SUBROUTINE F02SYF(N,D,E,NCOLB,B,LDB,NROWY,Y,LDY,NCOLZ,Z,LDZ,WORK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 14B REVISED. IER-839 (MAR 1990).
C  1. Purpose
C     =======
C
C  F02SYF  reduces  an  n by n  real bidiagonal matrix  to diagonal form
C  by  means  of  orthogonal transformations.  The  transformations  may
C  optionally be applied to given real matrices.
C
C  2. Description
C     ===========
C
C  The n by n bidiagonal matrix A is factorized as
C
C     A = Q*S*P',
C
C  where  Q and P  are  n by n orthogonal matrices and  S  is an  n by n
C  diagonal matrix  with  non-negative diagonal elements.  This  is  the
C  singular value decomposition  of the matrix A.
C
C  The diagonal elements of  S are the singular values of  A and will be
C  arranged in descending order. The columns of Q and P are the left and
C  right singular vectors of A respectively.
C
C  Optionally the matrices  C and/or Y and/or Z  given by
C
C     C = Q'*B,   W = Y*Q,   X = P'*Z,
C
C  where  B  is an  n by ncolb  matrix, Y  is an  nrowy by n  matrix and
C  Z  is an  n by ncolz  matrix, can also be returned.
C
C  The factorization is obtained by the  Golub-Reinsch version of the QR
C  algorithm.
C
C  3. Parameters
C     ==========
C
C  N      - INTEGER.
C
C           On entry, N specifies the order of the matrix  A.  N must be
C           at least  zero.  When  N = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  D      - REAL array of DIMENSION at least ( n ).
C
C           On entry, D must contain the diagonal elements of the n by n
C           bidiagonal matrix A such that
C
C              d( i ) = a( i, i ), i = 1, 2, ..., n.
C
C           On exit, D contains the  n singular values of  A arranged in
C           descending order of magnitude so that
C
C              d( 1 ) .ge. d( 2 ) .ge. ... .ge. d( n ) .ge. 0.
C
C  E      - REAL array of DIMENSION at least ( max( 1, n - 1 ) ).
C
C           On  entry,  E  must  contain  the  ( n - 1 )  super-diagonal
C           elements of the bidiagonal matrix A, with
C
C              e( i ) = a( i, i + 1 ), i = 1, 2, ..., n - 1.
C
C           E  is  used  as  internal  workspace  by  F02SYF  and so the
C           elements are changed on exit.
C
C  NCOLB  - INTEGER.
C
C           On entry,  NCOLB  must specify the  number of columns of the
C           matrix  B  and must be at least  zero.  When  NCOLB = 0  the
C           array B is not referenced.
C
C           Unchanged on exit.
C
C  B      - REAL             array of DIMENSION ( LDB, ncolb ).
C
C           Before entry with  NCOLB .gt. 0, the leading n by ncolb part
C           of the array B must contain the matrix to be transformed and
C           on exit  B  is overwritten by the  n by ncolb  matrix  Q'*B.
C
C           When  NCOLB = 0  the array  B  is not referenced.
C
C  LDB    - INTEGER.
C
C           On  entry,  LDB  must specify  the leading dimension  of the
C           array  B  as declared  in the  calling  (sub) program.  When
C           NCOLB .gt. 0  then LDB must be at least N.
C
C           Unchanged on exit.
C
C  NROWY  - INTEGER.
C
C           On entry,  NROWY  must specify  the  number of  rows  of the
C           matrix  Y  and must be at least  zero.  When  NROWY = 0  the
C           array Y is not referenced.
C
C           Unchanged on exit.
C
C  Y      - REAL             array of DIMENSION ( LDY, n ).
C
C           Before entry with  NROWY .gt. 0, the leading nrowy by n part
C           of the array Y must contain the matrix to be transformed and
C           on  exit  Y  is overwritten  by the  nrow by n  matrix  Y*Q.
C
C           When  NROWY = 0  the array  Y  is not referenced.
C
C  LDY    - INTEGER.
C
C           On  entry,  LDY  must specify  the leading dimension  of the
C           array  Y  as declared  in the  calling  (sub) program.  When
C           NROWY .gt. 0  then LDY must be at least NROWY.
C
C           Unchanged on exit.
C
C  NCOLZ  - INTEGER.
C
C           On entry,  NCOLZ  must specify the  number of columns of the
C           matrix  Z  and must be at least  zero.  When  NCOLZ = 0  the
C           array Z is not referenced.
C
C           Unchanged on exit.
C
C  Z      - REAL             array of DIMENSION ( LDZ, ncolz ).
C
C           Before entry with  NCOLZ .gt. 0, the leading n by ncolz part
C           of the array Z must contain the matrix to be transformed and
C           on exit  Z  is overwritten by the  n by ncolz  matrix  P'*Z.
C
C           When  NCOLZ = 0  the array  Z  is not referenced.
C
C  LDZ    - INTEGER.
C
C           On  entry,  LDZ  must specify  the leading dimension  of the
C           array  Z  as declared  in the  calling  (sub) program.  When
C           NCOLZ .gt. 0  then LDZ must be at least N.
C
C           Unchanged on exit.
C
C  WORK   - REAL array of DIMENSION at least ( max( 1, lwork ) ),  where
C           lwork must be at least zero when  ncolb = nrowy = ncolz = 0,
C           lwork must be at least   2*( n - 1 )   when  either  ncolb =
C           nrowy = 0  and  nrowz is positive, or  nrowz = 0  and one or
C           both of  ncolb and  nrowy  are positive, and  lwork  must be
C           at least  4*( n - 1 )  when one  or both of  ncolb and nrowy
C           are positive and  ncolz  is positive.
C
C           The array  WORK  is used as  internal workspace  by  F02SYF.
C           On exit,  WORK( 1 )  contains the total number of iterations
C           taken by the QR algorithm.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1  to specify noisy soft failure or noisy hard failure or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will be set to a  non-zero value  indicating either  that an
C           input parameter  has been  incorrectly set,  or that the  QR
C           algorithm  is not  converging.  See  the  next  section  for
C           further details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        N     .lt. 0
C        NCOLB .lt. 0
C        LDB   .lt. N      and  NCOLB .gt. 0
C        NROWY .lt. 0
C        LDY   .lt. NROWY  and  NROWY .gt. 0
C        NCOLZ .lt. 0
C        LDZ   .lt. N      and  NCOLZ .gt. 0
C
C  IFAIL .gt. 0
C
C     The  QR algorithm  has failed to converge in  50*N iterations.  In
C     this   case  d( 1 ), d( 2 ), ..., d( IFAIL )  may  not  have  been
C     found correctly  and the remaining singular values  may not be the
C     smallest.  The matrix  A will nevertheless have been factorized as
C     A = Q*E*P',  where  E is a bidiagonal matrix with  d( 1 ), d( 2 ),
C     ..., d( n )  as  the  diagonal elements  and  e( 1 ), e( 2 ), ...,
C     e( n - 1 )  as the super-diagonal elements.
C
C     This failure is not likely to occur.
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further Information
C     ===================
C
C  This routine  can be  used in conjunction  with routines  F02SWF  and
C  F02SXF   to  find  the  singular  value  decomposition  of  an  upper
C  triangular matrix.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 4-September-1987.
C     Sven Hammarling and Jeremy Du Croz, Nag Central Office.
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02SYF')
      DOUBLE PRECISION  QTR, ZERO
      PARAMETER         (QTR=0.25D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDB, LDY, LDZ, N, NCOLB, NCOLZ, NROWY
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LDB,*), D(*), E(*), WORK(*), Y(LDY,*),
     *                  Z(LDZ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAX, CS, DMAX, EKM2, SN, TEMP
      INTEGER           I, IERR, ITER, IW1, IW2, IW3, J, K, L, MAXIT, P
      LOGICAL           FORCE, WANTB, WANTY, WANTZ
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, F02XUS, F02XUT, F02XUU, F02XUV, F02XUW,
     *                  F06FGF, F06QKF, F06QXF, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      WORK(1) = ZERO
      WANTB = NCOLB .GT. 0
      WANTY = NROWY .GT. 0
      WANTZ = NCOLZ .GT. 0
      IERR = 0
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (NCOLB.LT.0) CALL P01ABY(NCOLB,'NCOLB',IFAIL,IERR,SRNAME)
      IF ((WANTB) .AND. (LDB.LT.N)) CALL P01ABY(LDB,'LDB',IFAIL,IERR,
     *    SRNAME)
      IF (NROWY.LT.0) CALL P01ABY(NROWY,'NROWY',IFAIL,IERR,SRNAME)
      IF ((WANTY) .AND. (LDY.LT.NROWY)) CALL P01ABY(LDY,'LDY',IFAIL,
     *    IERR,SRNAME)
      IF (NCOLZ.LT.0) CALL P01ABY(NCOLZ,'NCOLZ',IFAIL,IERR,SRNAME)
      IF ((WANTZ) .AND. (LDZ.LT.N)) CALL P01ABY(LDZ,'LDZ',IFAIL,IERR,
     *    SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
C
C     Find the size  of the  element  of  largest absolute value of  A.
C
      AMAX = ABS(D(1))
      DO 20 I = 2, N
         AMAX = MAX(AMAX,ABS(D(I)),ABS(E(I-1)))
   20 CONTINUE
C
C     Scale the matrix A by AMAX.
C
      IF (AMAX.GT.ZERO) THEN
         CALL DSCAL(N,1/AMAX,D,1)
         CALL DSCAL(N-1,1/AMAX,E,1)
      END IF
C
C     Split up the workspace.
C
      IW1 = 1
      IW2 = 1
      IW3 = 1
      IF ((WANTB) .OR. (WANTY)) THEN
         IW1 = N
         IF (WANTZ) THEN
            IW2 = (N-1) + IW1
            IW3 = (N-1) + IW2
         END IF
      ELSE IF (WANTZ) THEN
         IW3 = N
      END IF
C
C     Start the  QR algorithm.  A singular value is found for each value
C     of k.
C
      MAXIT = 50*N
      ITER = 1
      K = N
C+    WHILE( ( K.GT.1 ).AND.( ITER.LE.MAXIT ) )LOOP
   40 IF ((K.GT.1) .AND. (ITER.LE.MAXIT)) THEN
C
C        Test  to see if there is a  split,  or if we can force a split.
C        On exit from  F02XUT,  force will be  true  if we can  force  a
C        split and  p .gt. 0  means that  e( p )  is negligible, or will
C        be negligible after forcing the split.
C
         CALL F02XUT(ZERO,K,D,E,FORCE,P)
         L = P + 1
         IF (FORCE) THEN
            IF (P.EQ.K) THEN
C
C              Force a zero singular value in d( k ).
C
               CALL F02XUS(K,D,E,WANTZ,WORK(IW2),WORK(IW3))
               IF (WANTZ) CALL F06QXF('Left','Bottom','Backward',N,
     *                                NCOLZ,1,K,WORK(IW2),WORK(IW3),Z,
     *                                LDZ)
            ELSE
C
C              Force a split.
C
               CALL F02XUW(P,K,D,E,(WANTB) .OR. (WANTY),WORK,WORK(IW1))
               IF (WANTB) CALL F06QXF('Left','Top','Forward',N,NCOLB,P,
     *                                K,WORK,WORK(IW1),B,LDB)
               IF (WANTY) CALL F06QXF('Right','Top','Forward',NROWY,N,P,
     *                                K,WORK,WORK(IW1),Y,LDY)
            END IF
         END IF
         IF (L.GE.K) THEN
C
C           We have converged to a singular value.
C
            K = K - 1
         ELSE
C
C           Perform a QR step. First determine the shift.
C
            IF (K.GT.(L+1)) THEN
               EKM2 = E(K-2)
            ELSE
               EKM2 = ZERO
            END IF
            CALL F02XUU(-1,D(L),E(L),D(K-1),D(K),EKM2,E(K-1),CS,SN)
            CALL F02XUV('Non-zero',L,K,D,E,CS,SN,(WANTB) .OR. (WANTY),
     *                  WORK,WORK(IW1),WANTZ,WORK(IW2),WORK(IW3))
            IF (WANTB) CALL F06QXF('Left','Variable','Forward',N,NCOLB,
     *                             L,K,WORK,WORK(IW1),B,LDB)
            IF (WANTY) CALL F06QXF('Right','Variable','Forward',NROWY,N,
     *                             L,K,WORK,WORK(IW1),Y,LDY)
            IF (WANTZ) CALL F06QXF('Left','Variable','Forward',N,NCOLZ,
     *                             L,K,WORK(IW2),WORK(IW3),Z,LDZ)
            ITER = ITER + 1
         END IF
         GO TO 40
      END IF
C+    END WHILE
C
C     Unscale the matrix.
C
      IF (AMAX.GT.ZERO) THEN
         CALL DSCAL(N,AMAX,D,1)
         CALL DSCAL(N-1,AMAX,E,1)
      END IF
C
C     Make the singular values non-negative.
C
      DO 60 I = K, N
         IF (D(I).LT.ZERO) THEN
            D(I) = -D(I)
            IF (WANTB) CALL F06FGF(NCOLB,B(I,1),LDB)
            IF (WANTY) CALL F06FGF(NROWY,Y(1,I),1)
         END IF
   60 CONTINUE
C
C     Sort the singular values into descending order.
C
      DO 80 J = 1, K - 1
         WORK(J) = J + QTR
   80 CONTINUE
      DO 120 J = K, N
         DMAX = D(J)
         L = J
         DO 100 I = J + 1, N
            IF (D(I).GT.DMAX) THEN
               DMAX = D(I)
               L = I
            END IF
  100    CONTINUE
         WORK(J) = L + QTR
         IF (L.GT.J) THEN
            TEMP = D(J)
            D(J) = D(L)
            D(L) = TEMP
         END IF
  120 CONTINUE
      IF (WANTB) CALL F06QKF('Left','Transpose',N,WORK,NCOLB,B,LDB)
      IF (WANTY) CALL F06QKF('Right','No Transpose',N,WORK,NROWY,Y,LDY)
      IF (WANTZ) CALL F06QKF('Left','Transpose',N,WORK,NCOLZ,Z,LDZ)
      WORK(1) = ITER
      IF (K.EQ.1) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      ELSE
         WRITE (REC,FMT=99998) K
         IFAIL = P01ABF(IFAIL,K,SRNAME,2,REC)
      END IF
      RETURN
C
C
C     End of F02SYF. ( SBIQR )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
99998 FORMAT ('    The QR algorithm has failed to converge.',/'    ',I6,
     *       ' singular values have NOT been found.')
      END
