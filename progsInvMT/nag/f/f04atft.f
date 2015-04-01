      SUBROUTINE F04ATF(A,IA,B,N,C,AA,IAA,WKS1,WKS2,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Accurate solution of a set of real linear equations
C     with one right side.
C     1st. April 1973
C
C     Rewritten to call F07ADG, a modified version of LAPACK routine
C     SGETRF/F07ADF; new IFAIL exit inserted for illegal input
C     parameters; error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04ATF')
C     .. Scalar Arguments ..
      INTEGER           IA, IAA, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), AA(IAA,*), B(*), C(*), WKS1(*), WKS2(*)
C     .. Local Scalars ..
      INTEGER           I, IERR, IFAIL1, INFO, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04AHF, F07ADG
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99998) IA, N
      ELSE IF (IAA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99997) IAA, N
      ELSE
C
C        Copy A to working array AA
C
         DO 40 I = 1, N
            DO 20 J = 1, N
               AA(I,J) = A(I,J)
   20       CONTINUE
   40    CONTINUE
C
         CALL F07ADG(N,N,AA,IAA,WKS1,INFO)
C
         IF (INFO.EQ.0) THEN
            IF (N.GT.0) THEN
C
               IFAIL1 = 1
               CALL F04AHF(N,1,A,IA,AA,IAA,WKS1,B,N,X02AJF(),C,N,WKS2,N,
     *                     I,IFAIL1)
C
               IF (IFAIL1.NE.0) THEN
                  IERR = 2
                  NREC = 1
                  WRITE (P01REC,FMT=99995)
               END IF
C
            END IF
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99996)
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99997 FORMAT (1X,'** On entry, IAA.lt.max(1,N): IAA =',I16,', N =',I16)
99996 FORMAT (1X,'** Matrix A is approximately singular.')
99995 FORMAT (1X,'** Matrix A is too ill-conditioned.')
      END
