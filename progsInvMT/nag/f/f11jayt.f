      SUBROUTINE F11JAY(N,NNZ,NNZC,MAXF,A,IROW,ICOL,LFILL,DTOL,DSCALE,
     *                  ALPHA,PSTRAT,IPIV,NPIVM,NNZR,LEVF,IDLEVF,IR,IC,
     *                  ISTR,ISTC,ISTLL,LLNNZF,LLNNZB,IWORK,IERROR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     Computes a relaxed incomplete Cholesky factorization:
C
C               A  =  M + R
C
C                            T T
C               M  =  P L D L P
C
C     of a real NxN symmetric sparse matrix A, of arbitrary sparsity
C     pattern, stored in symmetric coordinate storage (SCS) format. In
C     the above decomposition L is lower triangular with unit diagonal
C     elements, D is diagonal, P is a permutation matrix and R is a
C     remainder matrix.
C
C     The matrix M is returned in terms of the SCS representation of
C     the lower triangular matrix
C
C                          -1
C               C  =  L + D  - I .
C
C     The amount of fill in the decomposition may be controlled either
C     by specifying the maximum level of fill LFILL, or by supplying a
C     drop tolerance DTOL below which fill-in elements are omitted.
C
C     The diagonal elements of A may be scaled prior to factorization
C     in order to ensure that M is positive-definite. Fill elements
C     which are not included in the factorization are multiplied by the
C     relaxation parameter ALPHA, and subtracted from the diagonal
C     element of the row currently being eliminated. Thus ALPHA = 0
C     corresponds to the standard incomplete Cholesky decomposition,
C     while ALPHA = 1 gives the modified incomplete Cholesky (MIC)
C     factorization, for which the row sums of the matrix are preserved.
C
C     The argument PSTRAT defines the pivoting strategy to be used. The
C     available options are no pivoting, user-defined diagonal pivoting,
C     and diagonal pivoting based on the Markowitz startegy, aimed at
C     minimizing fill-in.
C
C     Arguments
C     =========
C
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C            N >= 1.
C
C     NNZ    (input) INTEGER
C            On entry, the number of non-zero elements in the lower
C            triangular part of A.
C            0 <= NNZ <= N*(N+1)/2.
C
C     NNZC   (output) INTEGER
C            On exit, the number of non-zero elements in the lower
C            triangular matrix C.
C
C     MAXF   (input) INTEGER
C            On entry, the maximum amount of fill allowed. NNZ+MAXF
C            should be an overestimate of NNZC.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (NNZ+MAXF)
C            On entry, the non-zero elements in the lower triangular
C            part of the matrix A, ordered by increasing row index
C            and by increasing column index within each row. Multiple
C            entries for the same row and column indices are not
C            allowed.
C            On exit, A is overwritten by the non-zero elements of the
C            lower triangular matrix C.
C
C     IROW   (input/output) INTEGER array, dimension (NNZ+MAXF)
C     ICOL   (input/output) INTEGER array, dimension (NNZ+MAXF)
C            On entry, the row and column indices corresponding to the
C            non-zero elements given in the array A.
C            IROW and ICOL must satisfy the following constraints:
C            1 <= IROW(i) <= N, and 1 <= ICOL(i) <= IROW(i), for
C            i = 1,2,...,NNZ.
C            IROW(i-1) < IROW(i), or
C            IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i), for
C            i = 2,3,...,NNZ.
C            On exit, the row and column indices corresponding to the
C            non-zero elements returned in the array A.
C
C     LFILL  (input) INTEGER
C            On entry, the required level of fill.
C
C     DTOL   (input) DOUBLE PRECISION
C            On entry, if LFILL >= 0 DTOL is not used. Otherwise DTOL
C            gives the drop tolerance of the decomposition. A fill
C            element a   in row i and column j will be dropped if
C                     ij
C
C                |a  | < DTOL*SQRT(| a   * a   |).
C                  ij                 ii    jj
C
C            DTOL >= 0.0
C
C     DSCALE (input) DOUBLE PRECISION
C            On entry, the diagonal scaling parameter. Diagonal elements
C            of A are scaled by the factor (1.0 + DSCALE) prior to the
C            factorization.
C
C     ALPHA (input) DOUBLE PRECISION
C            On entry, the relaxation parameter. Fill elements which
C            are not included in the factorization are multiplied by
C            ALPHA and subtracted from the diagonal element of the row
C            currently being eliminated. For certain mesh-based problems
C            choosing ALPHA appropriately can reduce the order of
C            magnitude of the condition number of the preconditioned
C            matrix as a function of the mesh steplength.
C            0.0 <= ALPHA <= 1.0.
C
C     PSTRAT (input) CHARACTER*1
C            On entry, specifies the pivoting strategy to be adopted.
C               PSTRAT = 'N' or 'n' => no pivoting.
C               PSTRAT = 'U' or 'u' => user defined pivoting specified
C                                      by the input value of IPIV.
C               PSTRAT = 'M' or 'm' => diagonal pivoting aimed at
C                                      minimizing fill-in, based on
C                                      the Markowitz strategy.
C            PSTRAT = 'N', 'n', 'U', 'u', 'M' or 'm'.
C
C     IPIV   (input/output) INTEGER array, dimension (N)
C            On entry, if PSTRAT = 'U' or 'u' then IPIV(i) must specify
C            the row index of the diagonal element used as a pivot at
C            elimination stage i. Otherwise IPIV need not be defined.
C            On exit, the pivot indices. If IPIV(i) = j then the
C            diagonal element in row j was used as the pivot at
C            elimination stage i.
C
C     NPIVM  (output) INTEGER
C            On exit, the number of pivots which were modified during
C            the factorization to ensure that M is positive-definite.
C            If A is an M-matrix then no pivot modifications will be
C            required, and NPIVM = 0. For other cases the quality of
C            the preconditioner will generally depend on the returned
C            value of NPIVM. If NPIVM is large the preconditioner may
C            not be satisfactory. In this case it may be advantageous
C            to call F11JAY again with an increased value of DSCALE.
C
C     NNZR   (workspace) INTEGER array, dimension (N)
C            Used to store the number of non-zero elements in each row.
C
C     LEVF   (workspace) INTEGER array, dimension (MAXF+1)
C            Used to store the level of fill of each fill-in element.
C
C     IDLEVF (input) INTEGER
C            On entry, the dimension of the array LEVF as declared in
C            the (sub)program from which F11JAY is called.
C            IDLEVF >= MAXF+1 if LFILL >= 0, IDLEVF >= 1 otherwise.
C
C     IR     (workspace) INTEGER array, dimension (N)
C            Used to store row indices of non-zero elements in a given
C            column.
C
C     IC     (workspace) INTEGER array, dimension (NNZ+MAXF)
C            Used to store a copy of ICOL when it is converted to
C            linked-list format.
C
C     ISTR   (workspace) INTEGER array, dimension (N+1)
C            Used to store linked-list row start addresses. Holds
C            SCS storage row start addresses on exit.
C
C     ISTC   (workspace) INTEGER array, dimension (N)
C            Used to store linked-list column start addresses. Holds
C            SCS storage diagonal addresses on exit.
C
C     ISTLL  (workspace) INTEGER array, dimension (N)
C            Used to store start addresses for linked lists of rows
C            with equal numbers of non-zeros.
C
C     LLNNZF (workspace) INTEGER array, dimension (N)
C            Used to store forward linked lists of rows with equal
C            numbers of non-zeros.
C
C     LLNNZB (workspace) INTEGER array, dimension (N)
C            Used to store backward linked lists of rows with equal
C            numbers of non-zeros.
C
C     IWORK  (workspace) INTEGER array, dimension (N)
C            General workspace.
C
C     IERROR (output) INTEGER
C            On exit, the following values may occur:
C               IERROR = 0 => no error detected.
C               IERROR = 1 => MAXF is too small.
C               IERROR = 2 => error in internal call to F11ZBF.
C
C     ==================================================================
C     .. Parameters ..
      DOUBLE PRECISION  GROWTH
      PARAMETER         (GROWTH=0.01D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, DSCALE, DTOL
      INTEGER           IDLEVF, IERROR, LFILL, MAXF, N, NNZ, NNZC, NPIVM
      CHARACTER         PSTRAT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NNZ+MAXF)
      INTEGER           IC(NNZ+MAXF), ICOL(NNZ+MAXF), IPIV(N), IR(N),
     *                  IROW(NNZ+MAXF), ISTC(N), ISTLL(N), ISTR(N+1),
     *                  IWORK(N), LEVF(IDLEVF), LLNNZB(N), LLNNZF(N),
     *                  NNZR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, FAC, FSUB, FVAL, PIVMAX, PIVOT, RPIV,
     *                  SDPROD, SUM
      INTEGER           COL, I, IAFT, IBEF, ICI, IERYCF, INDMAX, IRI,
     *                  IWC, J, K, KENDR, L, LA, LCAUSE, LEAST, LELIM,
     *                  LEV, LL, M, NEXT, NNZRC, NNZRR, NR, PIVIND,
     *                  PIVROW, ROW
      LOGICAL           ADFILL, MARK
      CHARACTER         DUP, ZERO
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F06DFF, F11JAW, F11JAX, F11ZBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
C     Initialization.
C
      IERROR = 0
      NNZC = NNZ
      NPIVM = 0
C
C     Scale diagonal elements and move to the start of A, IROW and ICOL.
C
      LA = NNZ + MAXF
      CALL F11JAW(N,NNZC,A,LA,IROW,ICOL,DSCALE,IWORK,IERROR)
      IF (IERROR.EQ.1) RETURN
C
C     Set the logical array MARK to indicate whether Markowitz pivoting
C     is to be used. If necessary create doubly linked-lists used in the
C     selection of the Markowitz pivot.
C
      MARK = PSTRAT .EQ. 'M' .OR. PSTRAT .EQ. 'm'
C
      IF (MARK) THEN
C
C         Initialize NNZR, the number of non-zero elements in each row.
C
         CALL F06DBF(N,0,NNZR,1)
 
         DO 20 I = 1, NNZC
            IRI = IROW(I)
            ICI = ICOL(I)
            NNZR(IRI) = NNZR(IRI) + 1
            IF (IRI.NE.ICI) NNZR(ICI) = NNZR(ICI) + 1
   20    CONTINUE
C
C         Create doubly linked lists of rows which have equal numbers
C         of non-zero elements. These linked lists will be used in the
C         selection of pivot rows by the Markowitz process. LLNNZF and
C         LLNNZB hold the forward and backward pointers respectively,
C         while ISTLL holds the starting address of each list. Thus
C         ISTLL(i) holds the first address in LLNNZF of the list of
C         rows with i non-zeros.
C
         CALL F06DFF(N,NNZR,1,LLNNZF,1)
         CALL F11JAX(N,N,LLNNZF,ISTLL,IWORK)
C
         DO 40 I = 1, N
            J = LLNNZF(I)
            IF (J.GT.0) THEN
               LLNNZB(J) = I
            ELSE
               K = ISTLL(-J)
               LLNNZB(K) = J
            END IF
   40    CONTINUE
C
      ELSE
C
C         Use NNZR to store the pivot indices at each elimination stage.
C
         DO 60 I = 1, N
            NNZR(I) = I
   60    CONTINUE
         IF (PSTRAT.EQ.'U' .OR. PSTRAT.EQ.'u') THEN
            CALL F06DFF(N,IPIV,1,NNZR,1)
         END IF
C
      END IF
C
C     Convert SCS storage in IROW to linked-list format.
C
      CALL F11JAX(NNZC,N,IROW,ISTR,IWORK)
C
C     Copy ICOL into IC.
C
      CALL F06DFF(NNZC,ICOL,1,IC,1)
C
C     Convert SCS storage in ICOL to linked-list format.
C
      CALL F11JAX(NNZC,N,ICOL,ISTC,IWORK)
C
C     Set up IWORK as a dummy array to be used in the factorization.
C     Initialize the array IPIV of pivot indices to zero. A zero
C     entry indicates that the given row has not yet been used as
C     a pivot row.
C
      CALL F06DBF(N,1,IWORK,1)
      CALL F06DBF(N,0,IPIV,1)
C
C     Incomplete Cholesky factorization.
C
      DO 260 I = 1, N
C
C         Select pivot row.
C
         IF (MARK) THEN
C
C             Markowitz strategy.
C
            LEAST = 0
   80       CONTINUE
            LEAST = LEAST + 1
            LL = ISTLL(LEAST)
            IF (LL.EQ.0) GO TO 80
            J = LL
            PIVMAX = ABS(A(J))
            INDMAX = J
  100       CONTINUE
            AA = ABS(A(J))
            IF (AA.GT.PIVMAX .OR. (AA.EQ.PIVMAX .AND. J.LT.INDMAX)) THEN
               PIVMAX = AA
               INDMAX = J
            END IF
            K = LLNNZF(J)
            IF (K.GT.0) THEN
               J = K
               GO TO 100
            END IF
            PIVROW = INDMAX
C
C             Remove pivot row from the doubly linked lists.
C
            IBEF = LLNNZB(PIVROW)
            IAFT = LLNNZF(PIVROW)
            IF (IBEF.GT.0) THEN
               LLNNZF(IBEF) = IAFT
            ELSE
               ISTLL(-IBEF) = MAX(IAFT,0)
            END IF
            IF (IAFT.GT.0) LLNNZB(IAFT) = IBEF
C
         ELSE
C
            PIVROW = NNZR(I)
C
         END IF
C
         PIVIND = PIVROW
C
C         Reset IPIV(PIVROW) to indicate that the row PIVROW has been
C         chosen as a pivot row, and to store the permutation matrix.
C
         IPIV(PIVROW) = I
C
C         Scan the upper triangular part of the pivot row in the active
C         submatrix. If there is a non-zero element with address J in
C         column K set IWORK(K) = -J to indicate its presence. Also find
C         the number NR of such non-zero elements, store their column
C         indices IR and compute the sum of absolute values of non-zero
C         elements in the pivot row.
C
         NR = 0
         SUM = 0.D0
         J = ISTC(PIVROW)
C
  120    CONTINUE
C
C         Find column index.
C
         COL = J
  140    CONTINUE
         COL = IROW(COL)
         IF (COL.GT.0) GO TO 140
         COL = -COL
C
         IF (IPIV(COL).EQ.0) THEN
            IWORK(COL) = -J
            NR = NR + 1
            IR(NR) = COL
            SUM = SUM + ABS(A(J))
         END IF
C
         J = ICOL(J)
         IF (J.GT.0) GO TO 120
C
C         Scan the strictly lower triangular part of the pivot row in
C         the active submatrix marking non-zero entries in the dummy
C         array IWORK. Also update SUM, NR and IR to include these lower
C         triangular non-zero elements.
C
         J = ISTR(PIVROW)
  160    CONTINUE
         COL = IC(J)
C
         IF (IPIV(COL).EQ.0) THEN
            IWORK(COL) = -J
            NR = NR + 1
            IR(NR) = COL
            SUM = SUM + ABS(A(J))
         END IF
C
         J = IROW(J)
         IF (J.GT.0) GO TO 160
C
C         Modify pivot to ensure positive-definiteness and calculate
C         reciprocal pivot.
C
         PIVOT = A(PIVIND)
         IF (PIVOT.LE.GROWTH*SUM) THEN
            A(PIVIND) = SUM
            IF (SUM.EQ.0.D0) A(PIVIND) = 1.D0
            NPIVM = NPIVM + 1
         END IF
C
         RPIV = 1.D0/A(PIVIND)
C
C         Scan over non-pivot rows where elimination is required.
C         Find the row index, and the level of fill of the element
C         to be eliminated. Adjust the Markowitz counter array if
C         necessary to account for the zero being created.
C
         DO 220 M = 1, NR
C
            ROW = IR(M)
C
C             Calculate element of L used in elimination process.
C
            J = -IWORK(ROW)
C
            FAC = RPIV*A(J)
C
C             Find fill level of element being eliminated.
C
            LELIM = 0
            IF (J.GT.NNZ .AND. LFILL.GE.0) LELIM = LEVF(J-NNZ)
C
            IF (MARK) THEN
C
C                 Reduce NNZR(ROW) to account for the element being
C                 eliminated.
C
               NNZR(ROW) = NNZR(ROW) - 1
C
            END IF
C
C             For each entry in the non-pivot row check the entry in the
C             corresponding position in the dummy array IWORK. If it is
C             negative then both pivot and non-pivot rows have entries.
C             No fill will be generated, but the element is altered by
C             the elimination process. Change the sign of the dummy
C             element to show that this element has been visited.
C
            K = ISTR(ROW)
  180       CONTINUE
            COL = IC(K)
            IWC = IWORK(COL)
            IF (IWC.LT.0 .AND. IPIV(COL).EQ.0) THEN
               A(K) = A(K) - FAC*A(-IWC)
               IWORK(COL) = -IWC
            END IF
            KENDR = K
            K = IROW(K)
            IF (K.GT.0) GO TO 180
C
C             Scan the pivot row once again. For each entry check the
C             value in the corresponding position in the dummy array
C             IWORK. If it is positive, then this element has been
C             visited before, both pivot and non-pivot rows have an
C             element in this column, and hence there is no fill.
C             Change the sign of the dummy element to show that this
C             element has been revisited.
C             If the dummy array has a negative element there is an
C             entry in this column of the pivot row, but not in the
C             same column of the non-pivot row. There will therefore
C             be fill-in to this column of the non-pivot row. Add the
C             fill-in value to the end of the array A, increment NNZC,
C             update the linked-list pointers stored in IROW and ICOL,
C             and store the column index of the fill-in element in IC.
C
            DO 200 L = 1, NR
C
               COL = IR(L)
               IWC = IWORK(COL)
C
               IF (IWC.GT.0) THEN
C
                  IWORK(COL) = -IWC
C
               ELSE
C
                  IF (ROW.GE.COL) THEN
C
C                         Calculate fill-in and check whether it should
C                         be kept.
C
                     K = -IWC
                     FVAL = -FAC*A(K)
C
                     LCAUSE = 0
                     IF (K.GT.NNZ .AND. LFILL.GE.0) LCAUSE = LEVF(K-NNZ)
                     LEV = MAX(LELIM,LCAUSE) + 1
C
                     ADFILL = LEV .LE. LFILL
                     IF (LFILL.LT.0) THEN
                        SDPROD = SQRT(ABS(A(ROW)*A(COL)))
                        ADFILL = ABS(FVAL) .GE. DTOL*SDPROD
                     END IF
C
                     IF (ADFILL) THEN
C
C                             Add fill-in element.
C
                        NNZC = NNZC + 1
                        IF (NNZC-NNZ.GT.MAXF) THEN
                           IERROR = 1
                           RETURN
                        END IF
                        A(NNZC) = FVAL
                        IF (LFILL.GE.0) LEVF(NNZC-NNZ) = LEV
C
C                             If necessary adjust doubly linked lists
C                             of rows with equal numbers of non-zeros
C                             to take account of the change made to
C                             NNZR(COL).
C
                        IF (MARK) THEN
                           NNZR(ROW) = NNZR(ROW) + 1
                           NNZR(COL) = NNZR(COL) + 1
                           IBEF = LLNNZB(COL)
                           IAFT = LLNNZF(COL)
                           IF (IBEF.GT.0) THEN
                              LLNNZF(IBEF) = IAFT
                           ELSE
                              ISTLL(-IBEF) = MAX(IAFT,0)
                           END IF
                           IF (IAFT.GT.0) LLNNZB(IAFT) = IBEF
C
                           NNZRC = NNZR(COL)
                           LL = ISTLL(NNZRC)
                           IF (LL.NE.0) THEN
                              LLNNZF(COL) = LL
                              LLNNZB(LL) = COL
                           ELSE
                              LLNNZF(COL) = -NNZRC
                           END IF
                           ISTLL(NNZRC) = COL
                           LLNNZB(COL) = -NNZRC
                        END IF
C
C                             Make a new end to the IROW linked
C                             list and a new start to the ICOL
C                             linked list.
C
                        IROW(NNZC) = IROW(KENDR)
                        IROW(KENDR) = NNZC
                        KENDR = NNZC
                        ICOL(NNZC) = ISTC(COL)
                        ISTC(COL) = NNZC
C
C                             Store the column index of the fill
C                             element.
C
                        IC(NNZC) = COL
C
                     ELSE
C
C                             Subtract correction from diagonal
C
                        FSUB = ALPHA*FVAL
                        A(ROW) = A(ROW) + FSUB
                        A(COL) = A(COL) + FSUB
C
                     END IF
                  END IF
               END IF
C
  200       CONTINUE
C
C             If necessary adjust doubly linked lists of rows with
C             equal numbers of non-zeros to take account of the
C             changes made to NNZR(ROW).
C
            IF (MARK) THEN
               IBEF = LLNNZB(ROW)
               IAFT = LLNNZF(ROW)
               IF (IBEF.GT.0) THEN
                  LLNNZF(IBEF) = IAFT
               ELSE
                  ISTLL(-IBEF) = MAX(IAFT,0)
               END IF
               IF (IAFT.GT.0) LLNNZB(IAFT) = IBEF
C
               NNZRR = NNZR(ROW)
               LL = ISTLL(NNZRR)
               IF (LL.NE.0) THEN
                  LLNNZF(ROW) = LL
                  LLNNZB(LL) = ROW
               ELSE
                  LLNNZF(ROW) = -NNZRR
               END IF
               ISTLL(NNZRR) = ROW
               LLNNZB(ROW) = -NNZRR
            END IF
C
  220    CONTINUE
C
C         Reset dummy array.
C
         DO 240 M = 1, NR
            IWORK(IR(M)) = 1
  240    CONTINUE
C
  260 CONTINUE
C
C     Permute rows and columns into lower triangular order.
C
      CALL F06DFF(N,ISTR,1,IWORK,1)
C
      DO 280 I = 1, N
         ISTR(IPIV(I)) = IWORK(I)
  280 CONTINUE
C
      CALL F06DFF(N,ISTC,1,IWORK,1)
C
      DO 300 I = 1, N
         ISTC(IPIV(I)) = IWORK(I)
  300 CONTINUE
C
C     Convert linked-list representation back to SCS.
C
      DO 340 I = 1, N
         J = ISTR(I)
  320    CONTINUE
         NEXT = IROW(J)
         IROW(J) = I
         J = NEXT
         IF (NEXT.GT.0) GO TO 320
  340 CONTINUE
C
      DO 380 I = 1, N
         J = ISTC(I)
  360    CONTINUE
         NEXT = ICOL(J)
         ICOL(J) = I
         J = NEXT
         IF (NEXT.GT.0) GO TO 360
  380 CONTINUE
C
C     Exchange indices of any upper triangular elements.
C
      DO 400 I = N + 1, NNZC
         IRI = IROW(I)
         ICI = ICOL(I)
         IF (IRI.LT.ICI) THEN
            M = ICI
            ICOL(I) = IRI
            IROW(I) = M
         END IF
  400 CONTINUE
C                                            -1
C     Replace diagonal elements by those of D .
C
      DO 420 I = 1, N
         A(I) = 1.D0/A(I)
         IWORK(ICOL(I)) = I
  420 CONTINUE
C
C     Scale lower triangular elements by the diagonal for that column.
C
      DO 440 I = N + 1, NNZC
         A(I) = A(I)*A(IWORK(ICOL(I)))
  440 CONTINUE
C
C     Replace IPIV by its inverse.
C
      CALL F06DFF(N,IPIV,1,IWORK,1)
      DO 460 I = 1, N
         IPIV(IWORK(I)) = I
  460 CONTINUE
C
C     Sort non-zero elements into row order, with column order within
C     each row.
C
      DUP = 'R'
      ZERO = 'R'
      IERYCF = 0
      CALL F11ZBF(N,NNZC,A,IROW,ICOL,DUP,ZERO,ISTR,IWORK,IERYCF)
      IF (IERYCF.NE.0) IERROR = 2
C
      RETURN
      END
