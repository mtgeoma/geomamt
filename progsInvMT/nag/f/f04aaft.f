      SUBROUTINE F04AAF(A,IA,B,IB,N,M,C,IC,WKSPCE,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Approximate solution of a set of real linear equations
C     with multiple right hand sides.
C     1st August 1971
C
C     Rewritten to call F07ADG and F07AEG, modified versions of LAPACK
C     routines SGETRF/F07ADF and SGETRS/F07AEF; new IFAIL exit inserted
C     for illegal input parameters; error messages inserted.
C     February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AAF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IC, IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), B(IB,*), C(IC,*), WKSPCE(*)
C     .. Local Scalars ..
      INTEGER           I, IERR, INFO, J, NREC
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
      ELSE IF (M.LT.0) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) M
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99997) IA, N
      ELSE IF (IB.LT.MAX(1,N)) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99996) IB, N
      ELSE IF (IC.LT.MAX(1,N)) THEN
         IERR = 2
         NREC = 1
         WRITE (P01REC,FMT=99995) IC, N
      ELSE
C
C        Copy B to array C
C
         DO 40 I = 1, M
            DO 20 J = 1, N
               C(J,I) = B(J,I)
   20       CONTINUE
   40    CONTINUE
C
         CALL F07ADG(N,N,A,IA,WKSPCE,INFO)
C
         IF (INFO.EQ.0) THEN
C
            CALL F07AEG('No transpose',N,M,A,IA,WKSPCE,C,IC,INFO)
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99994)
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, M.lt.0: M =',I16)
99997 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99996 FORMAT (1X,'** On entry, IB.lt.max(1,N): IB =',I16,', N =',I16)
99995 FORMAT (1X,'** On entry, IC.lt.max(1,N): IC =',I16,', N =',I16)
99994 FORMAT (1X,'** Matrix A is approximately singular.')
      END
