      SUBROUTINE F04ARF(A,IA,B,N,C,WKS,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Approximate solution of a set of real linear
C     equations with one right hand side.
C     1st. April 1973
C
C     Rewritten to call F07ADG and F07AEG, modified versions of LAPACK
C     routines SGETRF/F07ADF and SGETRS/F07AEF; new IFAIL exit inserted
C     for illegal input parameters; error messages inserted.
C     February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04ARF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), B(*), C(*), WKS(*)
C     .. Local Scalars ..
      INTEGER           I, IERR, INFO, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F07ADG, F07AEG
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) IA, N
      ELSE IF (N.GT.0) THEN
C
C        Copy B to array C
C
         DO 20 I = 1, N
            C(I) = B(I)
   20    CONTINUE
C
         CALL F07ADG(N,N,A,IA,WKS,INFO)
C
         IF (INFO.EQ.0) THEN
C
            CALL F07AEG('No transpose',N,1,A,IA,WKS,C,N,INFO)
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99997)
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99997 FORMAT (1X,'** Matrix A is approximately singular.')
      END
