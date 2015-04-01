      SUBROUTINE F01LHF(N,NBLOKS,BLKSTR,A,LENA,PIVOT,TOL,INDEX,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 13B REVISED. IER-661 (AUG 1988).
C     MARK 14C REVISED. IER-884 (NOV 1990).
C
C     Purpose:
C
C     F01LHF decomposes a real almost block diagonal matrix.
C
C
C     Specification:
C
C        SUBROUTINE F01LHF(N,NBLOKS,BLKSTR,A,LENA,PIVOT,
C       *                  TOL,INDEX,IFAIL)
C     C     INTEGER    N,NBLOKS,BLKSTR(3,NBLOKS),LENA,PIVOT(N),
C     C    *           INDEX,IFAIL
C     C     real       A(LENA),TOL
C
C
C     Description:
C
C     The routine factorizes a real almost block diagonal matrix, A, by
C     column elimination with alternate row and column pivoting such
C     that no 'fill-in' is produced. The code, which is derived from
C     ARCECO desribed in [1], uses Level 2 BLAS [2].
C     No three successive diagonal blocks may have columns in
C     common and therefore the almost block diagonal matrix must have
C     the form shown in the following diagram.
C
C         .------|
C         | .    |
C     1   |___.__|
C           |   .-----------|
C           |     .         |
C           |       .       |
C           |         .     |
C     2     |___________.___|
C                   |     .-------|
C                   |       .     |
C     3             |_________.___|
C                               .
C                                 .
C                                   .
C                                     .
C                                       .
C                                   |-----.-----------|
C                                   |       .         |
C                                   |         .       |
C                                   |           .     |
C     NBLOKS-1                      |_____________.___|
C                                         |         .---|
C                                         |           . |
C     NBLOKS                              |_____________|
C
C
C     References:
C
C     [1] FORTRAN Packages for Solving Certain Almost Block Diagonal
C      Linear Systems by Modified Alternate Row and Column Elimination,
C      J.C.Diaz, G.Fairweather, P.Keast,
C      ACM Transactions on Mathematical Software (1983) 9, 358-375
C
C     [2] An Extended Set of Fortran Basic Linear Algebra Subprograms
C      J.J.Dongarra, J.DuCroz, S.Hammarling, R.J.Hanson
C      NAG Technical Report (1987) TR3/87
C
C
C     Parameters:
C
C       N - INTEGER
C           On entry, N must specify the number of rows of the
C           matrix A.
C           N .gt. 0.
C           Unchanged on exit.
C
C     NBLOKS - INTEGER
C           On entry, NBLOKS must specify the total number of blocks
C           of the matrix A.
C           0 .lt. NBLOKS .le. N.
C           Unchanged on exit.
C
C     BLKSTR - INTEGER array of DIMENSION (3,NBLOKS)
C           Before entry, BLKSTR must contain information which
C           describes the block structure of A as follows:
C           BLKSTR(1,k) must contain the number of rows in the kth
C                       block, k=1,2,..,NBLOKS;
C           BLKSTR(2,k) must contain the number of columns in the kth
C                       block, k=1,2,..,NBLOKS;
C           BLKSTR(3,k) must contain the number of columns of overlap
C                       between the kth and (k+1)th blocks,
C                       k=1,2,..,NBLOKS-1.
C           BLKSTR(3,NBLOKS) need not be set.
C
C           The following conditions delimit the structure of A.
C           BLKSTR(1,k),BLKSTR(2,k) .gt. 0  k=1,2,..,NBLOKS
C           BLKSTR(3,k) .ge. 0 k=1,2,..,NBLOKS
C           (there must be at least one column and one row in each block
C           and a non-negative number of columns of overlap)
C
C           BLKSTR(3,k-1)+BLKSTR(3,k) .le.
C                                   BLKSTR(2,k)  k=2,3,..,NBLOKS-1
C           (the total number of columns in overlaps in each block must
C           not exceed the number of columns in that block)
C
C           BLKSTR(2,1) .ge. BLKSTR(1,1)
C           BLKSTR(2,1)+sum(BLKSTR(2,k)-BLKSTR(3,k-1);k=2,j) .ge.
C                          sum(BLKSTR(1,k);k=1,j)  j=2,3,..,NBLOKS-1
C           sum(BLKSTR(2,k)-BLKSTR(3,k);k=1,j) .le.
C                          sum(BLKSTR(1,k);k=1,j)  j=1,2,..,NBLOKS-1
C           (the index of the first column of the overlap between the
C           jth and (j+1)th blocks must be .le. the index of the last
C           row of the jth block, and the index of the last column of
C           overlap must be .ge. the index of the last row of the jth
C           block)
C
C           sum(BLKSTR(1,k);k=1,NBLOKS) .eq. N
C           BLKSTR(2,1)+sum(BLKSTR(2,k)-BLKSTR(3,k-1);k=2,NBLOKS) .eq. N
C           (both the number of rows and of columns of A must equal N)
C
C           Unchanged on exit.
C
C       A - real array of DIMENSION (LENA)
C           Before entry, A must contain the elements of the real almost
C           block diagonal matrix stored block by block and each block
C           stored column by column. The sizes of the blocks and the
C           overlaps are defined by the parameter BLKSTR.
C           If the first element of the kth block is in position (m,n)
C           of the matrix A, then the (i,j)th element of the matrix A,
C           if it is in the kth block, should be stored in location
C               A(nbasek + (j+1-n)*nrowsk + (i+1-m))
C           where
C               nbasek = sum(BLKSTR(1,l)*BLKSTR(2,l);l=1,k-1)
C           is the base address of the kth block, and
C               nrowsk = BLKSTR(1,k)
C           which is the number of rows of the kth block.
C           On succesful exit, A contains the decomposed form of the
C           matrix.
C
C     LENA - INTEGER
C           On entry, LENA must specify the dimension of array A as
C           declared in the (sub)program from which F01LHF is called.
C           LENA .ge. sum(BLKSTR(1,k)*BLKSTR(2,k);k=1,NBLOKS).
C           Unchanged on exit.
C
C     TOL - real
C           Before entry, TOL must specify a relative tolerance to be
C           used to indicate whether or not the matrix is singular. See
C           'Further Comments' on how TOL is used. If TOL is non-
C           positive then TOL is reset to 10*eps, where eps is the
C           relative machine precision (see NAG Library routine X02AJF).
C
C     INDEX - INTEGER
C           On succesful exit with IFAIL = 0, or on exit with
C           IFAIL = -1, INDEX contains the value 0. On exit with
C           IFAIL = 1 INDEX contains the value k, where k is the
C           first position on the diagonal of the matrix A where too
C           small a pivot or a zero pivot was detected.
C
C     PIVOT - INTEGER array of DIMENSION (N)
C           On exit, PIVOT details the interchanges.
C
C     IFAIL - INTEGER
C
C           Before entry, IFAIL must be set to 0, -1 or 1. For users
C           not familiar with this parameter (described in chapter P01),
C           the recommended value is 0.
C           Unless the routine detects an error (see Section 6), IFAIL
C           contains 0 on exit.
C
C
C     Error Indicators:
C
C     Errors detected by the routine:
C
C     If on entry IFAIL = 0 or -1, explanatory error messages are
C     output on the current error message unit (as defined by NAG
C     Library routine X04AAF).
C
C     IFAIL = -1  N .lt. 1 or
C                 NBLOKS .lt. 1 or
C                 N .lt. NBLOKS or
C                 LENA too small or
C                 illegal values detected in BLKSTR.
C
C     IFAIL = 1   Decomposition complete but a small or zero pivot has
C                 been detected.
C
C
C     Further Comments:
C
C     Singularity or near singularity in A is determined by the
C     parameter TOL. If the absolute value of any pivot is less than
C     TOL * the maximum absolute element of A then A is said to be
C     singular. The position on the diagonal of A of the first of any
C     such pivots is indicated by the parameter INDEX.
C     This routine may be followed by  the routine F04LHF, which is
C     designed to solve sets of linear equations AX = B. In general,
C     more accurate solutions can be obtained by routine SUBSTB, but
C     are more expensive in terms of arithmetic operations and require
C     more storage.
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       AUXILIARY ROUTINES
C       ------------------
C       F01LHY: TO PERFORM THE COLUMN ELIMINATIONS WITH ROW PIVOTING
C       F01LHX: TO PERFORM THE COLUMN ELIMINATIONS WITH COLUMN PIVOTING
C       F01LHZ: TO PERFORM CHECKS ON THE PARAMETER BLKSTR
C       IDAMAX: TO FIND INDEX OF LARGEST ABSOLUTE ELEMEMT OF A VECTOR
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01LHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, INDEX, LENA, N, NBLOKS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LENA)
      INTEGER           BLKSTR(3,NBLOKS), PIVOT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  MAXELT, TL
      INTEGER           APTR, APTR2, CPVPTR, IDEV, IERR, LENARQ, M,
     *                  MXINDX, NCLBLK, NCLBOT, NCLPIV, NCPOLD, NOVRLP,
     *                  NRPVP1, NRWBLK, NRWBOT, NRWPIV, PIVPTR
      LOGICAL           REPORT
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           IDAMAX, P01ABF
      EXTERNAL          X02AJF, IDAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01LHX, F01LHY, F01LHZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      INDEX = 0
      IF (TOL.LE.0.0D0) TOL = 10.0D0*X02AJF()
      IERR = 0
      REPORT = IFAIL .LE. 0
      CALL X04AAF(0,IDEV)
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PERFORM DATA CONSISTENCY CHECKS .
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL F01LHZ(BLKSTR,NBLOKS,N,LENARQ,SRNAME,REPORT,IERR)
      IF (IERR.LT.0) THEN
         IERR = 1
      ELSE IF (LENARQ.GT.LENA) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC,FMT=99999) LENA, LENARQ
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (IERR.GT.0) GO TO 40
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     FIND LARGEST ELEMENT OF A .
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      MXINDX = IDAMAX(LENARQ,A,1)
      MAXELT = ABS(A(MXINDX))
      TL = TOL*MAXELT
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INTERNAL PARAMETERS
C
C      PIVPTR: POINTS TO FIRST ELEMENT IN A BLOCK OF THE ARRAY PIVOT
C              WHERE ROW INTERCHANGES WILL BE RECORDED .
C      CPVPTR: POINTS TO FIRST ELEMENT IN A BLOCK OF THE ARRAY PIVOT
C              WHERE COLUMN INTERCHANGES WILL BE RECORDED .
C      APTR  : POINTS TO START OF CURRENT BLOCK .
C      APTR2 : POINTS TO START OF NEXT BLOCK OF A .
C      NRWBLK: NUMBER OF ROWS IN CURRENT BLOCK .
C      NCLBLK: NUMBER OF COLUMNS IN CURRENT BLOCK.
C      NOVRLP: NOS OF COLUMNS OVERLAPPED BY CURRENT AND
C              NEXT BLOCKS .
C      NCLPIV: NUMBER OF COLUMN ELIMINATIONS WITH COLUMN PIVOTING
C              IN OVERLAP BETWEEN CURRENT AND NEXT BLOCKS .
C      NRWPIV: NUMBER OF COLUMN ELIMINATIONS WITH ROW PIVOTING
C              BEFORE OVERLAP BETWEEN CURRENT AND NEXT BLOCKS .
C      NCPOLD: NUMBER OF COLUMN ELIMINATIONS WITH COLUMN PIVOTING
C              IN OVERLAP BETWEEN LAST AND CURRENT BLOCKS .
C      NCLBOT: NUMBER OF COLUMNS IN NEXT BLOCK .
C      NRWBOT: NUMBER OF ROWS IN NEXT BLOCK .
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      APTR2 = 1
      NRWBOT = BLKSTR(1,1)
      NCLBOT = BLKSTR(2,1)
      NCLPIV = 0
      PIVPTR = 1
      CPVPTR = 1
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     LOOP OVER FIRST NBLOKS-1 BLOCKS PERFORMING :
C       NRWPIV COLUMN ELIMINATIONS WITH ROW PIVOTING FOLLOWED BY
C       NCLPIV COLUMN ELIMINATIONS WITH COLUMN PIVOTING
C     AND ON THE LAST BLOCK PERFORM NRWBLK(NBLOKS) COLUMN ELIMINATIONS
C     WITH ROW PIVOTING .
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO 20 M = 1, NBLOKS - 1
         APTR = APTR2
         NRWBLK = NRWBOT
         NCLBLK = NCLBOT
         NOVRLP = BLKSTR(3,M)
         NRWPIV = NCLBLK - NOVRLP - NCLPIV
         NRPVP1 = NRWPIV + 1
C
C        CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        CALL F01LHY TO PERFORM NRWPIV COLUMN ELIMS WITH ROW PIVOTING .
C
C        CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         CALL F01LHY(A(APTR),NRWBLK,NCLBLK,PIVOT(PIVPTR),NRWPIV,NCLPIV,
     *               TL,INDEX,IERR)
         PIVPTR = PIVPTR + NRWBLK
         CPVPTR = CPVPTR + NCLBLK - NOVRLP
         NCPOLD = NCLPIV
         NCLPIV = NRWBLK - NRWPIV
         NRWBOT = BLKSTR(1,M+1)
         NCLBOT = BLKSTR(2,M+1)
         APTR2 = APTR + NRWBLK*NCLBLK
C
C        CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        CALL F01LHX TO PERFORM NCLPIV COLUMN ELIMS WITH COLUMN PIVOTING
C
C        CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         CALL F01LHX(A(APTR2),NRWBOT,NCLBOT,A(APTR),NRWBLK,NCLBLK,
     *               PIVOT(CPVPTR),NCLPIV,NRWPIV,NCPOLD,NOVRLP,TL,INDEX,
     *               IERR)
         IF (IERR.EQ.0) INDEX = INDEX + NRWBLK
   20 CONTINUE
      NRWPIV = NCLBOT - NCLPIV
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     FINALLY PERFORM NRWBLK COLUMN ELIMS ON THE LAST BLOCK .
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL F01LHY(A(APTR2),NRWBOT,NCLBOT,PIVOT(PIVPTR),NRWPIV,NCLPIV,TL,
     *            INDEX,IERR)
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF (IERR.GT.0) THEN
         IERR = 2
         IF (REPORT) THEN
            WRITE (REC,FMT=99998) INDEX
            CALL X04BAF(IDEV,REC(1))
         END IF
      ELSE
         INDEX = 0
      END IF
C
   40 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
C
      RETURN
C
99999 FORMAT (' ** F01LHF - LENA(=',I16,') .lt. length required(=',I16,
     *       ') **')
99998 FORMAT (' ** F01LHF - too small a pivot detected for diagonal po',
     *       'sition ',I6)
      END
