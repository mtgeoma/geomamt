      SUBROUTINE Y90RLF(N,NBLOCK,BLOCK,SEED)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ===============================================
C         *  Y90RLF :  Generate random block structure  *
C         ===============================================
C
C
C              ARGUMENTS
C              =========
C
C     N      :  Integer, input.
C               Order of matrix to be generated.
C
C     NBLOCK :  Integer, output.
C               Number of blocks along the diagonal.
C
C     BLOCK  :  Integer array of DIMENSION (3,*), output.
C               For the i-th block BLOCK contains:
C               BLOCK(1,i)  ==>  Number of rows in the block.
C               BLOCK(1,i)  ==>  Number of columns in the block.
C               BLOCK(3,i)  ==>  Number of columns of overlap with the
C
C     SEED   :  Integer array of DIMENSION (4), Input/output.
C               Seeds for the random number generator.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
C     .. Scalar Arguments ..
      INTEGER           N, NBLOCK
C     .. Array Arguments ..
      INTEGER           BLOCK(3,*), SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  X1, X2, X3, XTOT
      INTEGER           I, IZERO, J, J1, J2, K1, K2, NBLMAX
C     .. External Functions ..
      DOUBLE PRECISION  Y90TBF
      EXTERNAL          Y90TBF
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F06DFF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, MOD, NINT, DBLE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Generate a random block-diagonal matrix
C
C-----------------------------------------------------------------------
C
C     Initialize all blocks to size = 2
C
      NBLOCK = N/2
      IF (N.LE.6) THEN
         NBLMAX = N/2
      ELSE IF (N.LE.10) THEN
         NBLMAX = N/3
      ELSE IF (N.LE.20) THEN
         NBLMAX = N/4
      ELSE
         NBLMAX = 6
      END IF
      CALL F06DBF(NBLOCK,2,BLOCK(1,1),3)
      BLOCK(1,NBLOCK/2+1) = BLOCK(1,NBLOCK/2+1) + MOD(N,2)
C
C     Randomly move lines (rows and columns) between blocks making sure
C     that no blocks have size lesser than 2.
C
      DO 40 I = 1, MAX(15,N)
         J1 = MAX(1,MIN(NINT(DBLE(NBLOCK)*Y90TBF(1,SEED)+HALF),NBLOCK))
         J2 = MOD(J1,NBLOCK) + 1
         K1 = BLOCK(1,J1)
         K2 = BLOCK(1,J2) + K1
         IF (BLOCK(1,J2).GE.NBLMAX) THEN
            BLOCK(1,J1) = BLOCK(1,J1) + 1
            BLOCK(1,J2) = BLOCK(1,J2) - 1
         ELSE IF (BLOCK(1,J1).GE.NBLMAX) THEN
            BLOCK(1,J2) = BLOCK(1,J2) + 1
            BLOCK(1,J1) = BLOCK(1,J1) - 1
         ELSE
            J = NINT(DBLE(K2)*Y90TBF(1,SEED)+HALF)
            IZERO = 0
            IF (J.LT.K1) THEN
               IF (BLOCK(1,J1).LE.2) THEN
                  IZERO = J1
                  BLOCK(1,J2) = BLOCK(1,J2) + BLOCK(1,J1)
                  BLOCK(1,J1) = 0
               ELSE
                  BLOCK(1,J2) = BLOCK(1,J2) + 1
                  BLOCK(1,J1) = BLOCK(1,J1) - 1
               END IF
            ELSE
               IF (BLOCK(1,J2).LE.2) THEN
                  IZERO = J2
                  BLOCK(1,J1) = BLOCK(1,J1) + BLOCK(1,J2)
                  BLOCK(1,J2) = 0
               ELSE
                  BLOCK(1,J1) = BLOCK(1,J1) + 1
                  BLOCK(1,J2) = BLOCK(1,J2) - 1
               END IF
            END IF
            IF (IZERO.GT.0) THEN
               J1 = BLOCK(1,IZERO)
               DO 20 J = IZERO, NBLOCK - 1
                  BLOCK(1,J) = BLOCK(1,J+1)
   20          CONTINUE
               NBLOCK = NBLOCK - 1
            END IF
         END IF
   40 CONTINUE
C
C     Make sure that there are no blocks of size 1
C
      J = 0
   60 CONTINUE
      J = J + 1
      IF (BLOCK(1,J).EQ.1) THEN
         J1 = MOD(J,NBLOCK) + 1
         BLOCK(1,J) = 0
         BLOCK(1,J1) = BLOCK(1,J1) + 1
         DO 80 I = J, NBLOCK - 1
            BLOCK(1,I) = BLOCK(1,I+1)
   80    CONTINUE
         NBLOCK = NBLOCK - 1
      END IF
      IF (J.LE.NBLOCK) GO TO 60
C-----------------------------------------------------------------------
C
C     Generate random overlaps between diagonal blocks
C
C-----------------------------------------------------------------------
C
C     1. Initialize
C
      CALL F06DBF(NBLOCK,0,BLOCK(3,1),3)
      CALL F06DFF(NBLOCK,BLOCK(1,1),3,BLOCK(2,1),3)
C
C     2. Generate the overlaps, such that at least a column of overlap
C        is generated for each contiguous pair of blocks and
C        non-contiguous blocks do not overlap.
C
      DO 100 I = 1, NBLOCK
         IF (I.LE.1) THEN
            J = NINT(DBLE(BLOCK(1,1)-1)*Y90TBF(1,SEED)+HALF)
            J = MAX(1,MIN(J,BLOCK(1,1)-1))
            BLOCK(2,2) = BLOCK(2,2) + J
            BLOCK(3,1) = J
         ELSE IF (I.LT.NBLOCK) THEN
            X1 = Y90TBF(1,SEED)
            X2 = Y90TBF(1,SEED)
            X3 = Y90TBF(1,SEED)
            XTOT = X1 + X2 + X3
            J1 = MAX(1,NINT(DBLE(BLOCK(1,I)-1)*X1/XTOT+HALF))
            J2 = MAX(1,NINT(DBLE(BLOCK(1,I)-1)*X2/XTOT+HALF))
            IF (J1+J2.GT.BLOCK(1,I)) THEN
               IF (J2.GT.J1) THEN
                  J2 = J2 - 1
               ELSE
                  J1 = J1 - 1
               END IF
            END IF
            BLOCK(2,I-1) = BLOCK(2,I-1) + J1
            BLOCK(2,I+1) = BLOCK(2,I+1) + J2
            BLOCK(3,I-1) = BLOCK(3,I-1) + J1
            BLOCK(3,I) = J2
         ELSE
            J = NINT(DBLE(BLOCK(1,I)-1)*Y90TBF(1,SEED)+HALF)
            J = MAX(1,MIN(J,BLOCK(1,I)-1))
            BLOCK(2,NBLOCK-1) = BLOCK(2,NBLOCK-1) + J
            BLOCK(3,NBLOCK-1) = BLOCK(3,NBLOCK-1) + J
         END IF
  100 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90RLF
C
C-----------------------------------------------------------------------
      RETURN
      END
