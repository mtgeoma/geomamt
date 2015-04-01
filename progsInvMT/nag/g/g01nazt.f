      SUBROUTINE G01NAZ(M,MPRTN,LDM,MR,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C       CONSTRUCTS THE MR X M MATRIX "MPRTN" CONTAINING
C       ALL MR PARTITIONS OF THE INTEGER M
C
C     .. Parameters ..
      INTEGER           MAXMOM
      PARAMETER         (MAXMOM=24)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDM, M, MR
C     .. Array Arguments ..
      INTEGER           MPRTN(LDM,M)
C     .. Local Scalars ..
      INTEGER           I, II, J, JJ, K, L, M1, N1, N2, N3
C     .. Local Arrays ..
      INTEGER           IWORK(MAXMOM), NUM(MAXMOM)
C     .. Data statements ..
      DATA              NUM/1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77,
     *                  101, 135, 176, 231, 297, 385, 490, 627, 792,
     *                  1002, 1255, 1575/
C     .. Executable Statements ..
C
      IF (M.LT.1 .OR. M.GT.MAXMOM) THEN
         IFAIL = 1
      ELSE IF (NUM(M).GT.LDM) THEN
         IFAIL = 2
      ELSE IF (M.EQ.1) THEN
         IFAIL = 0
         MR = 1
         MPRTN(1,1) = 1
      ELSE
         IFAIL = 0
         N1 = 0
         N2 = 1
         N3 = 0
         MR = 1
         M1 = 1
         L = 0
         MPRTN(1,1) = 1
         DO 20 J = 2, M
            MPRTN(1,J) = 0
   20    CONTINUE
         DO 220 K = 2, M
            IF (N2.NE.0) THEN
               DO 60 I = 1, N2
                  MPRTN(MR+I,1) = 0
                  DO 40 J = 2, M
                     MPRTN(MR+I,J) = MPRTN(I+N1,J)
   40             CONTINUE
                  MPRTN(MR+I,2) = MPRTN(MR+I,2) + 1
   60          CONTINUE
            END IF
            IF (N3.NE.0) THEN
               L = 0
               DO 180 I = N1 + N2 + 1, MR
                  DO 160 J = 2, K - 1
                     IF (MPRTN(I,J).NE.0) THEN
                        DO 80 JJ = 1, M
                           IWORK(JJ) = MPRTN(I,JJ)
   80                   CONTINUE
                        IWORK(J) = IWORK(J) - 1
                        IWORK(J+1) = IWORK(J+1) + 1
                        IF (N2.NE.0 .OR. L.NE.0) THEN
                           DO 120 II = MR + 1, MR + N2 + L
                              DO 100 JJ = 1, M - 1
                                 IF (MPRTN(II,JJ).NE.IWORK(JJ))
     *                               GO TO 120
  100                         CONTINUE
                              IF (MPRTN(II,M).EQ.IWORK(M)) GO TO 160
  120                      CONTINUE
                        END IF
                        L = L + 1
                        DO 140 JJ = 1, M
                           MPRTN(MR+N2+L,JJ) = IWORK(JJ)
  140                   CONTINUE
                     END IF
  160             CONTINUE
  180          CONTINUE
            END IF
            DO 200 II = 1, MR
               MPRTN(II,1) = MPRTN(II,1) + 1
  200       CONTINUE
            N1 = M1
            N3 = N2 + L
            N2 = MR - N1
            M1 = MR
            MR = MR + N3
  220    CONTINUE
      END IF
      RETURN
      END
