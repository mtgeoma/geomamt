      SUBROUTINE F03ADF(A,IA,N,DETR,DETI,WKSPCE,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Determinant of complex matrix
C     1st August 1971
C
C     Rewritten to call F07ARG, a modified version of LAPACK routine
C     CGETRF/F07ARF; new IFAIL exit inserted for illegal input
C     parameters; error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03ADF')
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DETI, DETR
      INTEGER           IA, IFAIL, N
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*)
      DOUBLE PRECISION  WKSPCE(*)
C     .. Local Scalars ..
      COMPLEX*16        DET, T
      DOUBLE PRECISION  ABSDET, ABST, OFLOW, UFLOW
      INTEGER           I, IERR, INFO, L, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F07ARG
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, MAX, NINT, DBLE
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 4
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 4
         NREC = 1
         WRITE (P01REC,FMT=99998) IA, N
      ELSE
C
         CALL F07ARG(N,N,A,IA,WKSPCE,INFO)
C
         IF (INFO.EQ.0) THEN
C
C           Compute the determinant as the product of the diagonal
C           elements of the factor L, with a factor + or - 1
C           determined by the interchanges.
C
            UFLOW = X02AMF()
            OFLOW = ONE/UFLOW
            DET = ONE
            DO 20 I = 1, N
               L = NINT(WKSPCE(I))
               IF (L.NE.I) DET = -DET
               T = A(I,I)
               ABST = ABS(DBLE(T)) + ABS(DIMAG(T))
               ABSDET = ABS(DBLE(DET)) + ABS(DIMAG(DET))
               IF (ABST.GE.ONE) THEN
                  IF (ABSDET.GT.OFLOW/ABST) THEN
                     IERR = 2
                     NREC = 1
                     WRITE (P01REC,FMT=99997)
                     DET = DCMPLX(OFLOW,OFLOW)
                     GO TO 40
                  ELSE
                     DET = DET*T
                  END IF
               ELSE
                  IF (ABSDET.LT.UFLOW/ABST) THEN
                     IERR = 3
                     NREC = 1
                     WRITE (P01REC,FMT=99996)
                     DET = ZERO
                     GO TO 40
                  ELSE
                     DET = DET*T
                  END IF
               END IF
   20       CONTINUE
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99995)
            DET = ZERO
         END IF
C
   40    DETR = DBLE(DET)
         DETI = DIMAG(DET)
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99997 FORMAT (1X,'** The value of the determinant is too large to be s',
     *       'tored.')
99996 FORMAT (1X,'** The value of the determinant is too small to be s',
     *       'tored.')
99995 FORMAT (1X,'** Matrix A is approximately singular.')
      END
