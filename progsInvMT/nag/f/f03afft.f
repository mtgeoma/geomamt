      SUBROUTINE F03AFF(N,EPS,A,IA,D1,ID,P,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     UNSYMDET
C     The unsymmetric matrix, A, is stored in the N*N array A(I,J),
C     I=1,N, J=1,N. The decomposition A=LU, where L is a
C     lower triangular matrix and U is a unit upper triangular
C     matrix, is performed and overwritten on A, omitting the unit
C     diagonal of U. A record of any interchanges made to the rows
C     of A is kept in P(I), I=1,N, such that the I-th row and
C     the P(I)-th row were interchanged at the I-th step. The
C     determinant, D1 * 2.0**ID, of A is also computed. The
C     subroutine will fail if A, modified by the rounding errors, is
C     singular or almost singular.
C     1st December 1971
C
C     Rewritten to call F07ADG, a modified version of LAPACK routine
C     SGETRF/F07ADF; new IFAIL exit inserted for illegal input
C     parameters; error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03AFF')
      DOUBLE PRECISION  ONE, ZERO, SIXTEN, SXTNTH
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0,SIXTEN=16.0D+0,
     *                  SXTNTH=6.25D-2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D1, EPS
      INTEGER           IA, ID, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), P(*)
C     .. Local Scalars ..
      INTEGER           I, IERR, INFO, L, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F07ADG
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, NINT
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
      ELSE
C
         CALL F07ADG(N,N,A,IA,P,INFO)
C
         IF (INFO.EQ.0) THEN
C
C           Compute the determinant as the product of the diagonal
C           elements of the factor L, with a factor + or - 1
C           determined by the interchanges.
C           Store the determinant as D = D1*2.0**ID.
C
            D1 = ONE
            ID = 0
            DO 60 I = 1, N
               L = NINT(P(I))
               IF (L.NE.I) D1 = -D1
               D1 = D1*A(I,I)
   20          IF (ABS(D1).GE.ONE) THEN
                  D1 = D1*SXTNTH
                  ID = ID + 4
                  GO TO 20
               END IF
   40          IF (ABS(D1).LT.SXTNTH) THEN
                  D1 = D1*SIXTEN
                  ID = ID - 4
                  GO TO 40
               END IF
   60       CONTINUE
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99997)
            D1 = ZERO
            ID = 0
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
