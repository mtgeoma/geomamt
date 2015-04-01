      SUBROUTINE D02TZF(MXMESH,NMESH,MESH,IPMESH,ERMX,IERMX,IJERMX,
     *                  RWORK,IWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02TZF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERMX
      INTEGER           IERMX, IFAIL, IJERMX, MXMESH, NMESH
C     .. Array Arguments ..
      DOUBLE PRECISION  MESH(MXMESH), RWORK(*)
      INTEGER           IPMESH(MXMESH), IWORK(*)
C     .. Local Scalars ..
      INTEGER           I, IER, IMSH, IPMSH, MXFIXP, N, NEQ, NFIXP,
     *                  NOLD, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
      IF (MXMESH.NE.IWORK(5)+1) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,i6,a,a/a,i6,a)')
     *     ' ** The value of MXMESH = ', MXMESH, ', which is not that',
     *     ' supplied to', ' ** the setup routine = ', IWORK(5) + 1, '.'
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
     *     , ' ** Error estimate information is not available.'
      ELSE IF (IWORK(16).EQ.6) THEN
         IER = 2
         NREC = 2
         WRITE (REC,FMT='(a/a)')
     *  ' ** The solver routine did not satisfy the error requirements.'
     *     , ' ** Information has been supplied on the last mesh used.'
      END IF
      IF (IER.EQ.1) GO TO 40
C 
C  Extract problem size dependent info from IWORK
C 
      NEQ = IWORK(2)
      MXFIXP = IWORK(20)
      NFIXP = IWORK(8)
C      IMSH = 1 + NEQ + MXFIXP
      NOLD = IWORK(11)
      NMESH = NOLD + 1
      N = IWORK(4)
C      IF (N.NE.NOLD) IMSH = IMSH + N + 1
      IMSH = 1 + MXFIXP + NEQ
      IPMSH = 20 + 3*NEQ
C
C
      DO 20 I = 1, NMESH
         MESH(I) = RWORK(IMSH+I)
         IPMESH(I) = IWORK(IPMSH+I)
   20 CONTINUE
      ERMX = RWORK(1)
      IERMX = IWORK(18)
      IJERMX = IWORK(19)
C
   40 CONTINUE
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
      END
