      SUBROUTINE E01SBF(M,X,Y,F,TRIANG,GRADS,PX,PY,PF,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Takes a Thiessen triangulation of a set of points in the plane,
C     and estimated partial derivatives with respect to X and Y, as
C     returned by NAG Fortran Library routine E01SAF. Returns the
C     value of a C1 function F(X,Y) interpolating the points and their
C     partial derivatives, evaluated at the point (PX,PY).
C     If (PX,PY) lies outside the triangulation boundary, extrapolation
C     is performed. Interpolation is exact for quadratic data.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01SBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PF, PX, PY
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  F(M), GRADS(2,M), X(M), Y(M)
      INTEGER           TRIANG(7*M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DUM1, DUM2
      INTEGER           IER, IST, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01SBZ
C     .. Save statement ..
      SAVE              IST
C     .. Data statements ..
      DATA              IST/1/
C     .. Executable Statements ..
      IER = 0
      NREC = 0
      CALL E01SBZ(M,PX,PY,X,Y,F,TRIANG(1),TRIANG(6*M+1),GRADS,IST,0,PF,
     *            DUM1,DUM2,IER)
      IF (IER.EQ.1) THEN
         NREC = 1
         WRITE (REC,FMT=99999) M
      ELSE IF (IER.EQ.2) THEN
         NREC = 2
         WRITE (REC,FMT=99998)
      ELSE IF (IER.EQ.3) THEN
         NREC = 2
         WRITE (REC,FMT=99997) PX, PY
      END IF
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, M .lt. 3: M =',I16,'.')
99998 FORMAT (1X,'** On entry, TRIANG does not contain a valid data po',
     *       'int triangulation;',/4X,'TRIANG may have been corrupted ',
     *       'since the call to E01SAF.')
99997 FORMAT (1X,'** Warning - the evaluation point (',1P,D12.4,',',
     *       D12.4,') lies outside the',/4X,'triangulation boundary. T',
     *       'he returned value was computed by extrapolation.')
      END
