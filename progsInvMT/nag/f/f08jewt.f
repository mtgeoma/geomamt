      SUBROUTINE F08JEW(ID,N,D,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  DLASRT sorts the numbers in D in increasing order
C  (if ID = 'I') or in decreasing order (if ID = 'D' ).
C
C  Use Quick Sort, reverting to Insertion sort on arrays of
C  size <= 20. Dimension of STACK limits N to about 2**32.
C
C  Arguments
C  =========
C
C  ID      (input) CHARACTER*1
C          = 'I': sort D in increasing order;
C          = 'D': sort D in decreasing order.
C
C  N       (input) INTEGER
C          The length of the array D.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, the array to be sorted.
C          On exit, D has been sorted into increasing order
C          (D(1) <= ... <= D(N) ) or into decreasing order
C          (D(1) >= ... >= D(N) ), depending on ID.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      INTEGER           SELECT
      PARAMETER         (SELECT=20)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
      CHARACTER         ID
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, D2, D3, DMNMX, TMP
      INTEGER           DIR, ENDD, I, J, START, STKPNT
C     .. Local Arrays ..
      INTEGER           STACK(2,32)
C     .. Executable Statements ..
C
C     Test the input paramters.
C
      INFO = 0
      DIR = -1
      IF (ID.EQ.'D' .OR. ID.EQ.'d') THEN
         DIR = 0
      ELSE IF (ID.EQ.'I' .OR. ID.EQ.'i') THEN
         DIR = 1
      END IF
      IF (DIR.EQ.-1) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      END IF
      IF (INFO.NE.0) THEN
C         CALL XERBLA('DLASRT',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.LE.1) RETURN
C
      STKPNT = 1
      STACK(1,1) = 1
      STACK(2,1) = N
   20 CONTINUE
      START = STACK(1,STKPNT)
      ENDD = STACK(2,STKPNT)
      STKPNT = STKPNT - 1
      IF (ENDD-START.LE.SELECT .AND. ENDD-START.GT.0) THEN
C
C        Do Insertion sort on D( START:ENDD )
C
         IF (DIR.EQ.0) THEN
C
C           Sort into decreasing order
C
            DO 60 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF (D(J).GT.D(J-1)) THEN
                     DMNMX = D(J)
                     D(J) = D(J-1)
                     D(J-1) = DMNMX
                  ELSE
                     GO TO 60
                  END IF
   40          CONTINUE
   60       CONTINUE
C
         ELSE
C
C           Sort into increasing order
C
            DO 100 I = START + 1, ENDD
               DO 80 J = I, START + 1, -1
                  IF (D(J).LT.D(J-1)) THEN
                     DMNMX = D(J)
                     D(J) = D(J-1)
                     D(J-1) = DMNMX
                  ELSE
                     GO TO 100
                  END IF
   80          CONTINUE
  100       CONTINUE
C
         END IF
C
      ELSE IF (ENDD-START.GT.SELECT) THEN
C
C        Partition D( START:ENDD ) and stack parts, largest one first
C
C        Choose partition entry as median of 3
C
         D1 = D(START)
         D2 = D(ENDD)
         I = (START+ENDD)/2
         D3 = D(I)
         IF (D1.LT.D2) THEN
            IF (D3.LT.D1) THEN
               DMNMX = D1
            ELSE IF (D3.LT.D2) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF (D3.LT.D2) THEN
               DMNMX = D2
            ELSE IF (D3.LT.D1) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
C
         IF (DIR.EQ.0) THEN
C
C           Sort into decreasing order
C
            I = START - 1
            J = ENDD + 1
  120       CONTINUE
  140       CONTINUE
            J = J - 1
            IF (D(J).LT.DMNMX) GO TO 140
  160       CONTINUE
            I = I + 1
            IF (D(I).GT.DMNMX) GO TO 160
            IF (I.LT.J) THEN
               TMP = D(I)
               D(I) = D(J)
               D(J) = TMP
               GO TO 120
            END IF
            IF (J-START.GT.ENDD-J-1) THEN
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = START
               STACK(2,STKPNT) = J
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = J + 1
               STACK(2,STKPNT) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = J + 1
               STACK(2,STKPNT) = ENDD
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = START
               STACK(2,STKPNT) = J
            END IF
         ELSE
C
C           Sort into increasing order
C
            I = START - 1
            J = ENDD + 1
  180       CONTINUE
  200       CONTINUE
            J = J - 1
            IF (D(J).GT.DMNMX) GO TO 200
  220       CONTINUE
            I = I + 1
            IF (D(I).LT.DMNMX) GO TO 220
            IF (I.LT.J) THEN
               TMP = D(I)
               D(I) = D(J)
               D(J) = TMP
               GO TO 180
            END IF
            IF (J-START.GT.ENDD-J-1) THEN
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = START
               STACK(2,STKPNT) = J
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = J + 1
               STACK(2,STKPNT) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = J + 1
               STACK(2,STKPNT) = ENDD
               STKPNT = STKPNT + 1
               STACK(1,STKPNT) = START
               STACK(2,STKPNT) = J
            END IF
         END IF
      END IF
      IF (STKPNT.GT.0) GO TO 20
      RETURN
C
C     End of F08JEW (DLASRT)
C
      END
