      SUBROUTINE G03ECZ(N,ILC,IUC,CD,IORD,DORD,IND,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes information for producing dendrogram
C
C     ILC(i) contains the cluster that joins IUC(i) at the
C     ith join. ILC(i) gt IUC(i)
C
C     IORD is the order for the dendrogram and object IORD(i)
C     joins cluster with object IORD(i-1) at distance DORD(i)
C
C     IND contains information on which part of the ILC is to be
C     searched
C
C     IERROR returns 1 if algorithm fails
C
C     .. Scalar Arguments ..
      INTEGER           IERROR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CD(N-1), DORD(N)
      INTEGER           ILC(N-1), IND(N), IORD(N), IUC(N-1)
C     .. Local Scalars ..
      INTEGER           I, ILINK, J, L
C     .. Executable Statements ..
C
      IERROR = 0
      DO 20 I = 1, N
         IND(I) = 1
   20 CONTINUE
      L = 1
      ILINK = 1
      IORD(1) = 1
      DORD(N) = CD(N-1)
      DO 140 I = 2, N
   40    CONTINUE
         DO 60 J = IND(L), N - 1
            IF (ILC(J).EQ.ILINK) THEN
               GO TO 100
            END IF
   60    CONTINUE
         IF (L.EQ.1) THEN
            IERROR = 1
            GO TO 160
         END IF
         IND(L) = N
   80    L = L - 1
         IF (IND(L).LT.N) THEN
            ILINK = IORD(L)
            GO TO 40
         ELSE IF (L.GT.1) THEN
            GO TO 80
         ELSE
            IERROR = 1
            GO TO 160
         END IF
  100    CONTINUE
         IND(L) = J + 1
         ILINK = IUC(J)
         IORD(I) = ILINK
         DORD(I-1) = CD(J)
  120    L = L + 1
         IF (L.GT.N) THEN
            IERROR = 1
            GO TO 160
         ELSE IF (IND(L).EQ.N) THEN
            GO TO 120
         END IF
  140 CONTINUE
  160 CONTINUE
      RETURN
      END
