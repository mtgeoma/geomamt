      SUBROUTINE F04LHF(TRANS,N,NBLOKS,BLKSTR,A,LENA,PIVOT,B,LDB,IR,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 13B REVISED. IER-665 (AUG 1988).
C     MARK 5 REVISED. IER-930 (APR 1991).
C
C     Purpose:
C
C     F04LHF calculates the approximate solution of a set of real linear
C     equations with multiple right hand sides, AX = B or A(trans)X = B,
C     where A is an almost block diagonal matrix which has been
C     decomposed using F01LHF.
C
C
C     Specification:
C
C        SUBROUTINE F04LHF(TRANS,N,NBLOKS,BLKSTR,A,LENA,PIVOT,
C       *                  B,LDB,IR,IFAIL)
C     C     INTEGER    N,NBLOKS,BLKSTR(3,NBLOKS),LENA,PIVOT(N),
C     C    *           LDB,IR,IFAIL
C     C     real       A(LENA),B(LDB,IR)
C     C     CHARACTER*1 TRANS
C
C
C     Description:
C
C     The routine solves a set of real linear equations AX = B or
C     A(trans)X = B, where
C     A is almost block diagonal. A must first be decomposed by
C     F01LHF. F04LHFthen computes X by forward and backward
C     substitution over the blocks. The code utilizes Level 3 BLAS [2].
C
C
C     References:
C
C     [1] FORTRAN Packages for Solving Certain Almost Block Diagonal
C      Linear Systems by Modified Alternate Row and Column Elimination,
C      J.C.Diaz, G.Fairweather, P.Keast,
C      ACM Transactions on Mathematical Software (1983) 9, 358-375
C
C     [2] A Proposal for a Set of Level 3 Basic Linear Algebra
C      Subprograms
C      J.J.Dongarra, J.DuCroz, I.Duff, S.Hammarling
C      NAG Technical Report (1987) In preparation
C
C
C     Parameters:
C
C     TRANS - CHARACTER*1.
C           On entry, TRANS must specify what solutions are required as
C           follows:
C                TRANS = 'N' or 'n'    AX = B
C                TRANS = 'T' or 't'    A(trans)X=B
C           Unchanged on exit.
C
C       N - INTEGER
C           On entry, N must specify the number of linear equations
C           defined by the matrix A. This must be the same parameter
C           N supplied to F01LHF when decomposing A.
C           N .gt. 0.
C           Unchanged on exit.
C
C     NBLOKS - INTEGER
C           On entry, NBLOKS must specify the total number of blocks
C           of the matrix A. This must be the same parameter NBLOKS
C           supplied to F01LHF.
C           0 .lt. NBLOKS .le. N.
C           Unchanged on exit.
C
C     BLKSTR - INTEGER array of DIMENSION (3,NBLOKS)
C           Before entry, BLKSTR must contain information which
C           describes the block structure of A. It must be unchanged
C           since the last call to F01LHF.
C           Unchanged on exit.
C
C       A - real array of DIMENSION (LENA)
C           Before entry, A must contain the elements in the
C           decomposition of A as output by F01LHF. It must be
C           unchanged since the last call to F01LHF.
C           Unchanged on exit.
C
C     LENA - INTEGER
C           On entry, LENA specifies the dimension of array A as
C           declared in the (sub)program from which F04LHF is called.
C           This must be the same parameter LENA supplied to F01LHF.
C           LENA .ge. sum(BLKSTR(1,k)*BLKSTR(2,k);k=1,NBLOKS)
C           Unchanged on exit.
C
C     PIVOT - INTEGER array of DIMENSION (N)
C           Before entry, PIVOT must contain the pivot information
C           about the decomposition as output from F01LHF. It must be
C           unchanged since the last call to F01LHF.
C           Unchanged on exit.
C
C       B - real array of DIMENSION (LDB,s) where s.ge.IR .
C           Before entry, B must contain the elements of the IR right
C           hand sides stored as columns.
C           On succesful exit, B will contain the IR solution vectors.
C
C     LDB - INTEGER
C           On entry, LDB must specify the first dimension of array B
C           as declared in the (sub)program from which F04LHF is called.
C           LDB .ge. N.
C           Unchanged on exit.
C
C      IR - INTEGER
C           On entry, IR must specify the number of right hand sides.
C           IR .gt. 0.
C           Unchanged on exit.
C
C     IFAIL - INTEGER
C           Before entry, IFAIL must be set to 0, -1 or 1. For users
C           not familiar with this parameter (described in chapter P01),
C           the recommended value is 0.
C           Unless the routine detects an error (see section 6), IFAIL
C           contains 0 on exit.
C
C
C     Error Indicators:
C
C     If on entry IFAIL = 0 or -1, explanatory error messages are output
C     on the current error message unit (as defined by the NAG Library
C     routine X04AAF).
C
C     IFAIL = 1  N .lt. 1 or
C                NBLOKS .lt. 1 or
C                IR .lt. 1 or
C                LDB .lt. N or
C                N .lt. NBLOKS or
C                LENA too small or
C                illegal values detected in BLKSTR or
C                illegal value for TRANS.
C
C
C     Further Comments:
C
C     None.
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       AUXILIARY ROUTINES
C       ------------------
C       SSWAP : TO INTERCHANGE TWO VECTORS
C       SGEMM : TO PERFORM MATRIX X MATRIX TYPE OPERATIONS
C       STRSM : TO PERFORM TRIANGULAR SOLVES
C       F04LHX: TO PERFORM TRAPEZOIDAL SOLVES
C       F01LHZ: TO PERFORM CHECKS ON THE PARAMETER BLKSTR
C
C     SGEMM AND STRSM ARE EXISTING LEVEL 3 BLAS
C     (THEY CALL SGEMV AND STRSV RESPECTIVELY)
C     F04LHX/STXSM IS A PROPOSED LEVEL 3 BLAS
C     (IT CALLS F04LHW/STXSV, A PROPOSED LEVEL 2 BLAS)
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     NOTE THAT THERE IS SPECIAL CODE TO HANDLE THE CASE WHEN IR = 1 .
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04LHF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IR, LDB, LENA, N, NBLOKS
      CHARACTER*1       TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LENA), B(LDB,IR)
      INTEGER           BLKSTR(3,NBLOKS), PIVOT(N)
C     .. Local Scalars ..
      INTEGER           IDEV, IERR, LENARQ
      LOGICAL           REPORT
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01LHZ, F04LHY, F04LHZ, X04AAF, X04BAF
C     .. Executable Statements ..
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INTERNAL PARAMETERS
C      NRWBLK: NUMBER OF ROWS IN A BLOCK .
C      NCLBLK: NUMBER OF COLUMNS IN A BLOCK .
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IERR = 0
      REPORT = IFAIL .LE. 0
      CALL X04AAF(0,IDEV)
      IF ((TRANS.NE.'N') .AND. (TRANS.NE.'t') .AND. (TRANS.NE.'T')
     *    .AND. (TRANS.NE.'n')) THEN
         IERR = -1
         IF (REPORT) THEN
            WRITE (REC,FMT=99996) TRANS
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (IR.LT.1) THEN
         IERR = -1
         IF (REPORT) THEN
            WRITE (REC,FMT=99999) IR
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      CALL F01LHZ(BLKSTR,NBLOKS,N,LENARQ,SRNAME,REPORT,IERR)
      IF (IERR.NE.-2) THEN
         IF (LDB.LT.N) THEN
            IERR = -1
            IF (REPORT) THEN
               WRITE (REC,FMT=99998) LDB, N
               CALL X04BAF(IDEV,REC(1))
            END IF
         END IF
      END IF
      IF (IERR.EQ.0) THEN
         IF (LENARQ.GT.LENA) THEN
            IERR = -1
            IF (REPORT) THEN
               WRITE (REC,FMT=99997) LENA, LENARQ
               CALL X04BAF(IDEV,REC(1))
            END IF
         END IF
      END IF
      IF (IERR.NE.0) THEN
         IERR = 1
         GO TO 20
      END IF
C
      IF (TRANS.EQ.'N' .OR. TRANS.EQ.'n') THEN
         CALL F04LHZ(N,NBLOKS,BLKSTR,A,LENA,PIVOT,B,LDB,IR)
      ELSE
         CALL F04LHY(N,NBLOKS,BLKSTR,A,LENA,PIVOT,B,LDB,IR)
      END IF
C
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
C
      RETURN
C
C
99999 FORMAT (' ** F04LHF - IR (=',I6,') .LT. 1 **')
99998 FORMAT (' ** F04LHF - LDB (=',I6,') .LT. N (=',I6,') **')
99997 FORMAT (' ** F04LHF - LENA (=',I6,') .LT. LENGTH REQUIRED (=',I6,
     *       ') **')
99996 FORMAT (' ** F04LHF - TRANS (=',A,') .NE. ''N'',''n'',''T'' OR ',
     *       '''t'' **')
      END
