      SUBROUTINE D02TYF(X,Y,NEQ,MMAX,RWORK,IWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C*****************************************************************
C     purpose
C           set up a standard call to  approx  to evaluate the
C           approximate solution  z = z( u(x) )  at a point x
C           (it has been computed by a call to  colnew ).
C           the parameters needed for  approx  are retrieved
C           from the work arrays  iwork  and  rwork .
C*****************************************************************
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02TYF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           IFAIL, MMAX, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(*), Y(NEQ,0:MMAX-1)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      INTEGER           I, ICOEF, IDMZ, IER, IMSH, ISPARE, IZ, J, K, KD,
     *                  MOMAX, MSTAR, MXFIXP, N, NMAX, NOLD, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  A(28), DUMMY(1)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           D02TKS, P01ABF
      EXTERNAL          D02TKS, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02TKT
C     .. Executable Statements ..
      IER = 0
      NREC = 0
C 
C  Check for valid call
C 
      IF (NEQ.NE.IWORK(2)) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,i6,a/a,i6,a)') ' ** The value of NEQ is ',
     *     NEQ, ' which is not that used in', ' ** the solver,',
     *     IWORK(2), '.'
      ELSE IF (IWORK(16).LT.1 .OR. IWORK(16).GT.6) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(a)')
     *     ' ** The solver routine does not appear to have been called.'
      ELSE IF (IWORK(16).GT.1 .AND. IWORK(16).LT.5) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a/a)')
     * ' ** The solver routine did not produce any results suitable for'
     *     , ' ** interpolation.'
      ELSE IF (IWORK(16).EQ.5) THEN
         IER = 2
         NREC = 3
         WRITE (REC,FMT='(a,a/a/a)')
     *     ' ** The solver routine did not converge to a suitable ',
     *     'solution.',
     *     ' ** An intermediate solution which converged has been used.'
     *     ,
     *   ' ** Interpolated values should be treated with great caution.'
      ELSE IF (IWORK(16).EQ.6) THEN
         IER = 2
         NREC = 2
         WRITE (REC,FMT='(a/a)')
     *  ' ** The solver routine did not satisfy the error requirements.'
     *     , ' ** Interpolated values should be treated with caution.'
      END IF
      IF (IER.EQ.1) GO TO 100
C 
C  Extract problem size dependent info from IWORK
C 
      K = IWORK(3)
      NMAX = IWORK(5)
      MSTAR = IWORK(7)
      MOMAX = IWORK(6)
      N = IWORK(4)
      NOLD = IWORK(11)
      MXFIXP = IWORK(20)
      IF (MOMAX.NE.MMAX) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,i1,a/a,i1,a)') ' ** The value of MMAX = ',
     *     MMAX, ' does not equal that specified',
     *     ' ** in the setup routine = ', MOMAX, '.'
         GO TO 100
      END IF
C 
C  Initialize Y
C 
      DO 40 I = 1, NEQ
         DO 20 J = 0, MMAX - 1
            Y(I,J) = 0.0D0
   20    CONTINUE
   40 CONTINUE
C 
C  Construct workspace pointers
C 
      KD = K*NEQ
      IMSH = 1 + NEQ + MXFIXP + 1
      IZ = IMSH + NMAX + 1
      IDMZ = IZ + MSTAR*(NMAX+1)
      ICOEF = IDMZ + KD*NMAX
      ISPARE = ICOEF + 49 + 2*(NMAX+1)
C 
C  Find interval containing the interpolation point
C 
      I = D02TKS(X,RWORK(IMSH),NOLD)
C 
C  Check for valid call
C 
      IF (I.LT.0) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,e11.4,a/a,e11.4,a,e11.4,a)')
     *     ' ** The value of X is ', X, ' which is outside the range',
     *     ' ** specified for the solver: ', RWORK(IMSH), ',',
     *     RWORK(IMSH+NOLD), '.'
         GO TO 100
      END IF
C 
C  Calculate interpolant values in 1-d array
C 
      CALL D02TKT(X,RWORK(ISPARE),A,RWORK(ICOEF),RWORK(IMSH+I-1),
     *            RWORK(IMSH+I),RWORK(IZ+(I-1)*MSTAR),RWORK(IDMZ+(I-1)
     *            *KD),K,NEQ,MMAX,IWORK(21),MSTAR,.TRUE.,.FALSE.,DUMMY)
C 
C  Extract interpolated values into Y
C 
      K = 0
      DO 80 I = 1, NEQ
         DO 60 J = 0, IWORK(20+I) - 1
            Y(I,J) = RWORK(ISPARE+K)
            K = K + 1
   60    CONTINUE
   80 CONTINUE
C 
  100 CONTINUE
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C 
      RETURN
      END
