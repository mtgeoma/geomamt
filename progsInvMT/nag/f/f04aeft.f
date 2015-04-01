      SUBROUTINE F04AEF(A,IA,B,IB,N,M,C,IC,WKSPCE,AA,IAA,BB,IBB,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Accurate solution of a set of real linear equations
C     with multiple right hand sides.
C     1st August 1971
C
C     Rewritten to call F07ADG, a modified version of LAPACK routine
C     SGETRF/F07ADF; new IFAIL exit inserted for illegal input
C     parameters; error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AEF')
C     .. Scalar Arguments ..
      INTEGER           IA, IAA, IB, IBB, IC, IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), AA(IAA,*), B(IB,*), BB(IBB,*), C(IC,*),
     *                  WKSPCE(*)
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
      ELSE IF (M.LT.0) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99998) M
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99997) IA, N
      ELSE IF (IB.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99996) IB, N
      ELSE IF (IC.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99995) IC, N
      ELSE IF (IAA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99994) IAA, N
      ELSE IF (IBB.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99993) IBB, N
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
         CALL F07ADG(N,N,AA,IAA,WKSPCE,INFO)
C
         IF (INFO.EQ.0) THEN
            IF (N.GT.0 .AND. M.GT.0) THEN
C
               IFAIL1 = 1
               CALL F04AHF(N,M,A,IA,AA,IAA,WKSPCE,B,IB,X02AJF(),C,IC,BB,
     *                     IBB,I,IFAIL1)
C
               IF (IFAIL1.NE.0) THEN
                  IERR = 2
                  NREC = 1
                  WRITE (P01REC,FMT=99991)
               END IF
C
            END IF
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99992)
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
99994 FORMAT (1X,'** On entry, IAA.lt.max(1,N): IAA =',I16,', N =',I16)
99993 FORMAT (1X,'** On entry, IBB.lt.max(1,N): IBB =',I16,', N =',I16)
99992 FORMAT (1X,'** Matrix A is approximately singular.')
99991 FORMAT (1X,'** Matrix A is too ill-conditioned.')
      END
