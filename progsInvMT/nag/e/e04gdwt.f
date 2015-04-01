      SUBROUTINE E04GDW(N,MAXRNK,S,IRANK)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-233 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04GDW PARTITIONS THE ORDERED ARRAY S, OF LENGTH N,
C     INTO TWO SECTIONS OF ROUGHLY EQUAL CONDITION NUMBER.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN AND
C     NICHOLAS I. M. GOULD.
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IRANK, MAXRNK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  S(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, S1, S2, SMXRNK
      INTEGER           I, IM1
C     .. Executable Statements ..
      IRANK = 1
      IF (MAXRNK.LE.2) RETURN
      S1 = S(1)
      S2 = S(2)
      SMXRNK = S(MAXRNK)
      A = 1.0D+0 + S2/SMXRNK
      IM1 = MAXRNK - 1
      DO 20 I = 2, IM1
         B = S1/S(I) + S(I+1)/SMXRNK
         IF (B.GE.A) GO TO 20
         IRANK = I
         A = B
   20 CONTINUE
      RETURN
C
C     END OF E04GDW   (PARTN)
C
      END
