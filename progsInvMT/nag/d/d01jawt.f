      SUBROUTINE D01JAW(IRAD2)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     D01JAW GENERATES, AT EACH CALL, A NEW COMBINATION
C     IX(1)...IX(N) SATISFYING
C         IX(1).GE.IX(2).GE. ... .GE.IX(N).GE.0,    (1)
C
C         IX(1)**2 + ... + IX(N)**2 = IRAD2,        (2)
C     AND
C         NOT ALL IX(I)  I=1...N  EVEN.             (3)
C
C     THE COMBINATIONS ARE GENERATED IN INVERSE-LEXICOGRAPHIC
C     ORDER (SEE COMMENTS IN THE ROUTINE D01JAX).
C
C     MAJOR VARIABLES
C     ---------------
C
C     SEE ALSO COMMENTS IN THE ROUTINE D01JAZ.
C
C     IX2    - IX2(I) = IX(I)**2   I=1...N.
C
C     INDEX  - ON ENTRY, INDEX IS THE LARGEST NUMBER SUCH THAT
C              IX(INDEX).GT.0.
C              EXCEPTIONS- A) IX(N).GT.0, INDEX = N-1.
C                          B) INDEX = 0, IX IS UNDEFINED. D01JAW IS
C                                      CALLED FOR THE FIRST TIME
C                                      WITH THE GIVEN VALUE OF
C                                      IRAD2.
C
C     IREST  - AT EACH STAGE OF THE COMPUTATION,
C              IREST = IX2(INDEX+1) + ... + IX(N).
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IRAD2
C     .. Scalars in Common ..
      INTEGER           INDEX, INOG, IREST, N, NMIN, NPLUS
C     .. Arrays in Common ..
      INTEGER           IX(4), IX2(4)
C     .. Local Scalars ..
      INTEGER           I, IND1, IVALUE, NEXT, NEXT2
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT, INT
C     .. Common blocks ..
      COMMON            /BD01JA/N, NMIN, NPLUS, IX, IX2, INDEX, INOG,
     *                  IREST
C     .. Executable Statements ..
      IF (INDEX.NE.0) GO TO 40
      IREST = IRAD2
      INDEX = 1
      DO 20 I = 2, N
         IX(I) = 0
   20 CONTINUE
      NEXT = IREST
      GO TO 100
C
C     IX CONTAINS A GOOD COMBINATION (SATISFYING THE CONDITIONS (1),
C     (2) AND (3)). GENERATE THE NEXT GOOD COMBINATION.
C
C     SET IX(INDEX) TO IX(INDEX)-1, UNLESS THIS WOULD PREVENT A NEW
C     GOOD COMBINATION BEING GENERATED.
C
   40 IREST = IREST + IX2(INDEX)
      NEXT = IX(INDEX) - 1
      NEXT2 = NEXT*NEXT
      IF ((NPLUS-INDEX)*NEXT2.GE.IREST) GO TO 60
      INDEX = INDEX - 1
      IF (INDEX.EQ.0) GO TO 220
      GO TO 40
C
C     SET IX(INDEX) AND COMPUTE NEW VALUE OF IX(I)  I.GT.INDEX.
C
   60 IX(INDEX) = NEXT
      IX2(INDEX) = NEXT2
      IREST = IREST - NEXT2
C
C     COMPUTE THE COMBINATION IX(1)...IX(N-1), SUCH THAT IX(1)..
C     ..IX(INDEX) IS LEFT UNCHANGED AND IX(1)**2+...+IX(N-1)**2
C     IS AS LARGE AS POSSIBLE BUT .LE.IRAD2.
C
   80 IF (INDEX.EQ.NMIN) GO TO 120
      INDEX = INDEX + 1
  100 IVALUE = INT(SQRT(DBLE(IREST)+0.1D0))
      IF (IVALUE.GT.NEXT) IVALUE = NEXT
      IX2(INDEX) = IVALUE*IVALUE
      IX(INDEX) = IVALUE
      IREST = IREST - IX2(INDEX)
      IF (IREST.EQ.0) GO TO 140
      NEXT = IVALUE
      GO TO 80
C
C     COMPUTE, IF POSSIBLE, THE VALUE OF IX(N) SUCH THAT
C     IX(1)...IX(N-1)IX(N) IS A GOOD COMBINATION, AND EXIT.
C
  120 IVALUE = INT(SQRT(DBLE(IREST)+0.1D0))
      IF (IREST.NE.IVALUE*IVALUE) GO TO 40
      IX(N) = IVALUE
      GO TO 180
  140 IND1 = INDEX + 1
      DO 160 I = IND1, N
         IX(I) = 0
  160 CONTINUE
  180 DO 200 I = 1, N
         IF (IX(I).NE.IX(I)/2*2) RETURN
  200 CONTINUE
      GO TO 40
  220 INOG = 0
      RETURN
      END
