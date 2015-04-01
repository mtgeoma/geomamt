      SUBROUTINE E04GDU(N,S,IRANK)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-233 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04GDU REPARTITIONS THE REAL ARRAY S SO THAT IRANK IS
C     REDUCED.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN,
C     BRIAN T. HINDE AND ENID M. LONG
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C
C     IF THE RANK OF THE ARRAY IS ALREADY 0, RETURN TO THE MAIN
C     PROGRAM.
C
C     .. Scalar Arguments ..
      INTEGER           IRANK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  S(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, EPSRNK, PROD
      INTEGER           IR1, K
C     .. Intrinsic Functions ..
      INTRINSIC         LOG10, INT
C     .. Executable Statements ..
      IF (IRANK.EQ.0) RETURN
C
C     IF THE RANK OF THE ARRAY IS 1, REDUCE THE RANK TO 0.
C
      IF (IRANK.GT.1) GO TO 20
      IRANK = 0
      RETURN
C
C     SET EPSRNK EQUAL TO THE POWER OF 10.0 LARGER THAN THE
C     IRANK TH ELEMENT OF S.
C
   20 IR1 = IRANK
      A = S(IR1)
      PROD = LOG10(A)
      K = INT(PROD) + 1
      EPSRNK = 1.0D+1**K
C
C     CHECK TO SEE IF ANY OTHER ELEMENTS OF S ARE OF THE SAME
C     ORDER OF MAGNITUDE AS THE (IRANK - 1) ST. IF SO MODIFY THE
C     RANK APPROPRIATELY.
C
   40 IR1 = IR1 - 1
      IRANK = IR1
      IF (IR1.EQ.0) RETURN
      A = S(IR1)
      IF (A.LT.EPSRNK) GO TO 40
      RETURN
C
C     END OF E04GDU   (RDRANK)
C
      END
