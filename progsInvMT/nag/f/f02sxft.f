      SUBROUTINE F02SXF(N,A,LDA,NCOLY,Y,LDY,WORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F02SXF returns either the matrix Z given by
C
C     Z = P'*Y,
C
C  where  Y is an  n by ncoly matrix and  P is the right-hand orthogonal
C  matrix  associated  with the  transformation of an  upper  triangular
C  matrix to bidiagonal form, or the matrix  P'  itself.
C
C  This routine must be preceded by a call to routine F02SWF.
C
C  2. Description
C     ===========
C
C  Routine  F02SWF factorizes the  n by n upper triangular matrix  R  as
C
C     R = Q*B*P',
C
C  where  Q and P  are  n by n orthogonal matrices and  B  is an  n by n
C  bidiagonal matrix and information on the matrix  P is returned in the
C  upper triangular part of the array A. Following a call to F02SWF this
C  routine may be used to find either the  n by ncoly  matrix  Z  or the
C  n by n orthogonal matrix P'.
C
C  3. Parameters
C     ==========
C
C  N      - INTEGER.
C
C           On entry, N specifies the order of the matrix  R.  N must be
C           at least  zero.  When  N = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n )
C
C           Before entry, the leading  N by N  upper triangular  part of
C           the array  A  must contain information on the matrix  P,  as
C           returned  by  routine  F02SWF.
C
C           On exit  with   NCOLY .gt. 0,  the  array  A  is  unchanged.
C           On exit  with   NCOLY = 0,  the leading  N by N  part of the
C           array  A  is overwritten by the matrix  P'.
C
C  LDA    - INTEGER.
C
C           On  entry,  LDA  must specify  the leading dimension  of the
C           array  A as declared in the calling (sub) program. LDA  must
C           be at least N.
C
C           Unchanged on exit.
C
C  NCOLY  - INTEGER.
C
C           On entry ,  NCOLY  must specify the number of columns of the
C           matrix  Y  and must be at least  zero. When  NCOLY = 0  then
C           the array  Y  is not referenced, but instead the array  A is
C           overwritten by the matrix  P'.
C
C           Unchanged on exit.
C
C  Y      - REAL             array of DIMENSION ( LDY, ncoly ).
C
C           Before entry with  NCOLY .gt. 0, the leading n by ncoly part
C           of the array Y must contain the matrix to be transformed and
C           on exit  Y  is overwritten by the  n by ncoly  matrix  P'*Y.
C
C           When  NCOLY = 0  then the array  Y is not referenced.
C
C
C  LDY    - INTEGER.
C
C           On  entry,  LDY  must specify  the leading dimension  of the
C           array  Y  as  declared  in the  calling (sub) program.  When
C           NCOLY .gt. 0  then LDY must be at least n.
C
C           Unchanged on exit.
C
C  WORK   - REAL             array    of      DIMENSION      at    least
C           max( 1, 2*( n - 1 ) ).
C
C           Used as internal workspace.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to  -1  indicating that an  input parameter has
C           been  incorrectly  set. See  the  next  section  for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        N     .lt. 0
C        LDA   .lt. N
C        NCOLY .lt. 0
C        NCOLY .gt. 0  and  LDY .lt. N
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C
C  Nag auxiliary linear algebra routine.
C
C  -- Written on 22-July-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02SXF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDY, N, NCOLY
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(*), Y(LDY,*)
C     .. Local Scalars ..
      INTEGER           IERR, J, K
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06BCF, F06FBF, F06QXF, P01ABY
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.N) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (NCOLY.LT.0) CALL P01ABY(NCOLY,'NCOLY',IFAIL,IERR,SRNAME)
      IF ((NCOLY.GT.0) .AND. (LDY.LT.N)) CALL P01ABY(LDY,'LDY',IFAIL,
     *    IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      IF (NCOLY.GT.0) THEN
C
C        Form Z as
C
C           Z = P( n - 2 )*...*P( 2 )*P( 1 )*Y.
C
         DO 40 K = 1, N - 2
C
C           Recover the rotations that put the zeros into the kth row of
C           R. The cosines and sines that define  P( k ) are returned in
C           work( 1 ), ..., work( n - 1 )       and      work( n ), ...,
C           work( 2*( n - 1 ) ).
C
            DO 20 J = K + 2, N
               CALL F06BCF(A(K,J),WORK(J-1),WORK(N+J-2))
   20       CONTINUE
C
C           Form P( k )*Y.
C
            CALL F06QXF('Left','Variable pivot','Backwards',N,NCOLY,K+1,
     *                  N,WORK,WORK(N),Y,LDY)
   40    CONTINUE
      ELSE
C
C        Form P' as
C
C           P' = I*( P( n - 2 )*...*P( 2 )*P( 1 ) ).
C
         IF (N.GT.1) THEN
            A(N,N) = ONE
            A(N-1,N) = ZERO
            A(N,N-1) = ZERO
            IF (N.GT.2) THEN
               DO 80 K = N - 2, 1, -1
                  A(K+1,K+1) = ONE
                  A(K,K+1) = ZERO
                  DO 60 J = K + 2, N
                     CALL F06BCF(-A(K,J),WORK(J-1),WORK(N+J-2))
                     A(K,J) = ZERO
   60             CONTINUE
                  CALL F06FBF(N-K,ZERO,A(K+1,K),1)
                  CALL F06QXF('Right','Variable pivot','Forward',N-K,
     *                        N-K,1,N-K,WORK(K+1),WORK(N+K),A(K+1,K+1),
     *                        LDA)
   80          CONTINUE
            END IF
         END IF
         A(1,1) = ONE
      END IF
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F02SXF. ( SBIAP  )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
