      SUBROUTINE F11JAW(N,NNZ,A,LA,IROW,ICOL,DSCALE,IWORK,IERROR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     Given a symmetric coordinate storage (SCS) representation of a
C     real sparse symmetric matrix A F11JAW scales the diagonal elements
C     and moves them to the first N storage locations. If there is no
C     diagonal element for a given row the routine adds a zero element
C     there.
C
C     Arguments
C     =========
C
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C
C     NNZ    (input/output) INTEGER
C            On entry, the number of non-zero elements in the lower
C            triangular part of A.
C            On exit, the number of non-zero elements in the lower
C            triangular part of A, including any zero diagonals added.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (LA)
C            On entry, the non-zero elements in the lower triangular
C            part of the matrix A, ordered by increasing row index
C            and by increasing column index within each row. Multiple
C            entries for the same row and column indices are not
C            allowed.
C            On exit, the reordered non-zero elements, including any
C            zero diagonal elements added.
C
C     LA     (input) INTEGER
C            On entry, the dimension of the arrays A, IROW and ICOL
C            as declared in the (sub)program from which F11JAW is
C            called.
C
C     IROW   (input/output) INTEGER array, dimension (LA)
C            On entry, the row indices corresponding to the non-zero
C            elements given in the array A.
C            On exit, the row indices corresponding to the non-zero
C            elements returned in the array A.
C
C     ICOL   (input/output) INTEGER array, dimension (LA)
C            On entry, the column indices corresponding to the non-zero
C            elements given in the array A.
C            On exit, the column indices corresponding to the non-zero
C            elements returned in the array A.
C
C     DSCALE (input) DOUBLE PRECISION
C            On entry, the diagonal scaling parameter. Diagonal elements
C            of A are scaled by the factor (1 + DSCALE).
C
C     IWORK  (workspace) INTEGER array, dimension (N)
C
C     IERROR (output) INTEGER
C            On exit, the following values may occur:
C               IERROR = 0 => no error detected.
C               IERROR = 1 => LA is too small.
C
C     ==================================================================
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DSCALE
      INTEGER           IERROR, LA, N, NNZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA)
      INTEGER           ICOL(LA), IROW(LA), IWORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACPY, DS
      INTEGER           I, ICI, IRI
C     .. External Subroutines ..
      EXTERNAL          F06DBF
C     .. Executable Statements ..
C
C     Add any missing diagonal entries.
C
      CALL F06DBF(N,0,IWORK,1)
C
      DO 20 I = 1, NNZ
         IRI = IROW(I)
         ICI = ICOL(I)
         IF (IRI.EQ.ICI) IWORK(IRI) = IWORK(IRI) + 1
   20 CONTINUE
C
      IERROR = 0
      DO 40 I = 1, N
         IF (IWORK(I).EQ.0) THEN
            NNZ = NNZ + 1
            IF (NNZ.GT.LA) THEN
               IERROR = 1
               RETURN
            END IF
            A(NNZ) = 0.D0
            IROW(NNZ) = I
            ICOL(NNZ) = I
         END IF
   40 CONTINUE
C
C     Move diagonal entries to the first N locations.
C
      I = 1
   60 CONTINUE
C
      IRI = IROW(I)
      ICI = ICOL(I)
      IF ((IRI.EQ.ICI .AND. IRI.NE.I) .AND. (IROW(IRI)
     *    .NE.IRI .OR. ICOL(IRI).NE.IRI)) THEN
         ACPY = A(I)
         A(I) = A(IRI)
         A(IRI) = ACPY
         IROW(I) = IROW(IRI)
         ICOL(I) = ICOL(IRI)
         IROW(IRI) = IRI
         ICOL(IRI) = IRI
      ELSE
         IF (I.EQ.NNZ) GO TO 80
         I = I + 1
      END IF
C
      GO TO 60
C
   80 CONTINUE
C
C     Scale diagonal entries.
C
      DS = 1.D0 + DSCALE
      DO 100 I = 1, N
         A(I) = DS*A(I)
  100 CONTINUE
C
      RETURN
      END
