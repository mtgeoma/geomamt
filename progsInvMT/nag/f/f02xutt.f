      SUBROUTINE F02XUT(TEST,N,D,E,FORCE,P)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  F02XUT  tests a bidiagonal matrix to see if there is a split, or if a
C  split can be forced.
C
C  The bidiagonal matrix must be supplied in the arrays D and E with the
C  diagonal  elements  in  d( 1 ), d( 2 ), ..., d( n )  and  the  super-
C  diagonal  elements  in  e( 1 ), e( 2 ), ..., e( n - 1 ).  The  arrays
C  are unaltered on exit.
C
C  The  bidiagonal matrix  is  searched  from the  bottom backwards  for
C  negligible elements  and  P  returns  the  row  in  which  the  first
C  negligible element is found. If no negligible element is found then
C  P = 0  is returned.
C
C  If  d( p ) is negligible then  FORCE  is set to true, indicating that
C  a  split  can be forced,  otherwise  FORCE  is returned as false.  If
C  FORCE  is  returned as  false  and  P  is  returned  as positive then
C  e( p ) is negligible.
C
C  When   test .le. 0.0   then  a  local test  is  used  and  d( p )  is
C  regarded as negligible if
C
C    abs( d( p ) ) .le. eps*( max( abs( e( p ) ), abs( e( p - 1 ) ) ) ),
C
C  where  e( n ) = e( 0 ) = 0.0   and   eps   is  the  relative  machine
C  precision as returned by routine  X02AJF, or if
C
C    abs( d( p ) ) .lt. flmin/eps,
C
C  where  flmin  is the  underflow  threshold  as  returned  by  routine
C  X02AKF, and similarly  e( p )  is regarded as negligible if
C
C    abs( e( p ) ) .le. eps*( max( abs( d( p + 1 ) ), abs( d( p ) ) ) ),
C
C  or if
C
C    abs( e( p ) ) .lt. flmin/eps.
C
C  When   test .gt. 0.0  then  a  global test  is  used  and  d( p )  is
C  regarded as negligible if
C
C    abs( d( p ) ) .lt. eps*test,
C
C  and similarly  e( p )  is regarded as negligible if
C
C    abs( e( p ) ) .lt. eps*test.
C
C  n  must be at least zero.
C
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 1-October-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TEST
      INTEGER           N, P
      LOGICAL           FORCE
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSD, ABSDI, ABSE, ABSEI, AMAX, DMAX, EMAX, EPS,
     *                  FLMIN, NEGL, SMALL
      INTEGER           I
      LOGICAL           FIRST
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AKF
      EXTERNAL          X02AJF, X02AKF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Save statement ..
      SAVE              FIRST, EPS, FLMIN, SMALL
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
      IF (FIRST) THEN
         FIRST = .FALSE.
         EPS = X02AJF()
         FLMIN = X02AKF()
         SMALL = FLMIN/EPS
      END IF
C
      FORCE = .FALSE.
      I = N
      IF (TEST.LE.ZERO) THEN
         IF (N.EQ.1) THEN
            IF (ABS(D(N)).LT.SMALL) THEN
               FORCE = .TRUE.
               GO TO 60
            END IF
         ELSE
            ABSD = ABS(D(N))
            ABSE = ABS(E(N-1))
            AMAX = MAX(ABSD,ABSE)
            IF ((ABSD.LE.EPS*AMAX) .OR. (AMAX.LT.SMALL)) THEN
               FORCE = .TRUE.
               GO TO 60
            END IF
            DO 20 I = N - 1, 2, -1
               ABSDI = ABS(D(I))
               DMAX = MAX(ABSDI,ABSD)
               AMAX = MAX(DMAX,ABSE)
               IF ((ABSE.LE.EPS*DMAX) .OR. (AMAX.LT.SMALL)) GO TO 60
               ABSEI = ABS(E(I-1))
               EMAX = MAX(ABSE,ABSEI)
               AMAX = MAX(EMAX,ABSDI)
               IF ((ABSDI.LE.EPS*EMAX) .OR. (AMAX.LT.SMALL)) THEN
                  FORCE = .TRUE.
                  GO TO 60
               END IF
               ABSD = ABSDI
               ABSE = ABSEI
   20       CONTINUE
            ABSDI = ABS(D(1))
            DMAX = MAX(ABSDI,ABSD)
            AMAX = MAX(DMAX,ABSE)
            IF ((ABSE.LE.EPS*DMAX) .OR. (AMAX.LT.SMALL)) GO TO 60
            AMAX = MAX(ABSE,ABSDI)
            IF ((ABSDI.LE.EPS*ABSE) .OR. (AMAX.LT.SMALL)) THEN
               FORCE = .TRUE.
               GO TO 60
            END IF
         END IF
      ELSE
         NEGL = EPS*TEST
         IF (ABS(D(N)).LT.NEGL) THEN
            FORCE = .TRUE.
            GO TO 60
         END IF
         DO 40 I = N - 1, 1, -1
            IF (ABS(E(I)).LT.NEGL) GO TO 60
            IF (ABS(D(I)).LT.NEGL) THEN
               FORCE = .TRUE.
               GO TO 60
            END IF
   40    CONTINUE
      END IF
      I = 0
   60 CONTINUE
      P = I
      RETURN
C
C     End of F02XUT. ( SBITST )
C
      END
