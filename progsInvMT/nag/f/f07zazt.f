      SUBROUTINE F07ZAZ(ISPEC,NAME,IVAL,RWFLAG)
C
C  Mark 15 Release. NAG Copyright 1991.
C  Mark 16 Revised. NAG Copyright 1993.
C  Mark 17 Revised. NAG Copyright 1995.
C
C  This version of F07ZAZ was generated
C  automatically by program GENZAZ.
C
C  -- NAG version of LAPACK auxiliary routine ILAENV
C
C  The parameter MAXIC is set to the largest value
C  that can be returned by F07ZAY. NCODES is the
C  number of Fortran Library routines that are
C  expected to call F07ZAZ. NSPECS is the number
C  of values that ISPEC can take.
C
C  Purpose
C  =======
C
C  F07ZAZ sets or returns problem-dependent
C  parameters for the local environment. See
C  ISPEC for a description of the parameters.
C
C  The problem-dependent parameters are contained
C  in the integer array IPARMS, and the value with
C  index ISPEC is set or copied to IVAL.
C
C  Arguments
C  =========
C
C  ISPEC (input) INTEGER
C     Specifies the parameter to be set or
C     returned by F07ZAZ.
C     = 1: the optimal blocksize; if this value
C          is 1, an unblocked algorithm will give
C          the best performance.
C     = 2: the minimum block size for which the
C          block routine should be used; if the
C          usable block size is less than this
C          value, an unblocked routine should be
C          used.
C     = 3: the crossover point (in a block
C          routine, for N less than this value,
C          an unblocked routine should be used).
C     = 4: the number of shifts, used in the
C          nonsymmetric eigenvalue routines.
C     = 5: the minimum column dimension for
C          blocking to be used; rectangular
C          blocks must have dimension at least
C          k by m, where k is given by
C          F07ZAZ(2,...) and m by F07ZAZ(5,...).
C     = 6: the crossover point for the SVD (when
C          reducing an m by n matrix to bidiagonal
C          form, if max(m,n)/min(m,n) exceeds this
C          value, a QR factorization is used first
C          to reduce the matrix to a triangular
C          form).
C     = 7: the number of processors.
C     = 8: the crossover point for the multishift
C          QR and QZ methods for nonsymmetric
C          eigenvalue problems.
C
C  NAME  (input) CHARACTER*(*)
C     The name of the calling subroutine.
C
C  IVAL  (input/output) INTEGER
C     the value of the parameter set or returned.
C
C  FLAG  (input) INTEGER
C     = 0: F07ZAZ returns in IVAL the value of
C          the parameter specified by ISPEC.
C     = 1: F07ZAZ sets the parameter specified
C          by ISPEC to the value in IVAL.
C
C  ==============================================
C
C     .. Parameters ..
      INTEGER           NSPECS, NCODES, MAXIC
      PARAMETER        (NSPECS=8,NCODES=  48,MAXIC= 124)
C     .. Scalar Arguments ..
      INTEGER           ISPEC, IVAL, RWFLAG
      CHARACTER*(*)     NAME
C     .. Local Scalars ..
      INTEGER           I, ICODE
C     .. Local Arrays ..
      INTEGER           IPARMS(NSPECS,NCODES),
     +                  POINT(0:MAXIC)
C     .. External Functions ..
      INTEGER           F07ZAY
      EXTERNAL          F07ZAY
C     .. Save statement ..
      SAVE              IPARMS
C     .. Data statements ..
      DATA              (IPARMS(I,  1),I=1,NSPECS)
     +                  /  16,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  2),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  3),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  4),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  5),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  6),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  7),I=1,NSPECS)
     +                  /  16,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  8),I=1,NSPECS)
     +                  /  16,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I,  9),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 10),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 11),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 12),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 13),I=1,NSPECS)
     +                  /  16,   5,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 14),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 15),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 16),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 17),I=1,NSPECS)
     +                  /   1,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 18),I=1,NSPECS)
     +                  /   8,   2, 224,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 19),I=1,NSPECS)
     +                  /   8,   2, 192,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 20),I=1,NSPECS)
     +                  /   8,   2,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 21),I=1,NSPECS)
     +                  /   8,   2, 320,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 22),I=1,NSPECS)
     +                  /   8,   2, 256,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 23),I=1,NSPECS)
     +                  /   8,   2,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 24),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 25),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 26),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 27),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 28),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 29),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 30),I=1,NSPECS)
     +                  /   8,   2, 192,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 31),I=1,NSPECS)
     +                  /   8,   2,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 32),I=1,NSPECS)
     +                  /   8,   2, 280,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 33),I=1,NSPECS)
     +                  /   8,   2,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 34),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 35),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 36),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 37),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 38),I=1,NSPECS)
     +                  /   8,   2, 448,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 39),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 40),I=1,NSPECS)
     +                  /   0,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 41),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 42),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 43),I=1,NSPECS)
     +                  /   8,   2, 192,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 44),I=1,NSPECS)
     +                  /   1,   1,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 45),I=1,NSPECS)
     +                  /   0,   0,   0,   3,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 46),I=1,NSPECS)
     +                  /   0,   0,   0,   2,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 47),I=1,NSPECS)
     +                  /  32,   0,   0,   0,   0,   0,   0,   0 /
      DATA              (IPARMS(I, 48),I=1,NSPECS)
     +                  /  24,   0,   0,   0,   0,   0,   0,   0 /
      DATA        POINT /
     +                   0,   1,   2,   3,   4,   5,   0,   6,
     +                   0,   7,   8,   9,  10,  11,   0,  12,
     +                   0,  13,   0,  14,   0,   0,   0,  15,
     +                   0,   0,  16,   0,  17,  18,  19,  20,
     +                  21,  22,  23,  24,  25,  26,  27,  28,
     +                  29,   0,  30,  31,  32,   0,  33,   0,
     +                  34,  35,  36,   0,  37,  38,   0,   0,
     +                   0,   0,   0,  39,   0,   0,   0,   0,
     +                   0,   0,   0,   0,   0,  40,   0,   0,
     +                   0,   0,   0,   0,   0,  41,   0,   0,
     +                   0,   0,   0,  42,   0,   0,   0,   0,
     +                   0,  43,   0,   0,   0,   0,   0,  44,
     +                   0,   0,   0,   0,   0,  45,   0,   0,
     +                   0,   0,   0,  46,   0,   0,   0,   0,
     +                   0,  47,   0,   0,   0,   0,   0,  48,
     +                   0,   0,   0,   0,   0
     +                  /
C     .. Executable Statements ..
C
C     Convert the NAG name to an integer code.
      ICODE = POINT(F07ZAY(NAME))
C
      IF (ISPEC.LT.1 .OR. ISPEC.GT.NSPECS) THEN
C        Invalid value for ISPEC
         IVAL = -1
      ELSE IF (ICODE.EQ.0) THEN
C        Invalid value for NAME
         IVAL = -2
      ELSE IF (RWFLAG.EQ.0) THEN
C        Read the value of a parameter
         IVAL = IPARMS(ISPEC,ICODE)
      ELSE
C        Set the value of a parameter
         IPARMS(ISPEC,ICODE) = IVAL
      END IF
C
      RETURN
C
C     End of F07ZAZ
C
      END
