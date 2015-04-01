      SUBROUTINE E04HEY(N,LH,LPH,NS,IGRADE,V,IV,P,PHESL,PHESD,PRHS,H,W)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     **************************************************************
C
C     E04HEY FORMS THE PROJECTIONS VTHV AND VTHP WHERE T DENOTES THE
C     TRANSPOSE. P IS A VECTOR AND THE SYMMETRIC MATRIX H IS STORED
C     BY ROWS AS A LOWER-TRIANGULAR MATRIX IN THE VECTOR H OF LENGTH
C     N*(N + 1)/2. ON OUTPUT, THE SYMMETRIC MATRIX VT * H * V IS
C     STORED IN THE FORM PHESL(TRANSPOSE) * PHESD * PHESL, WHERE
C     PHESL IS A LOWER TRIANGULAR MATRIX WITH UNIT DIAGONALS STORED
C     AS A STRICT LOWER TRIANGLE OF LENGTH LPH AND PHESD IS A
C     DIAGONAL MATRIX OF LENGTH NS.
C
C     NOTE THAT THE WORKSPACE ARRAY W MUST BE OF LENGTH AT LEAST N.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN AND
C     NICHOLAS I. M. GOULD.
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     Modified to call BLAS.
C     Peter Mayes, NAG Central Office, October 1987.
C
C     .. Scalar Arguments ..
      INTEGER           IGRADE, IV, LH, LPH, N, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LH), P(N), PHESD(NS), PHESL(LPH), PRHS(NS),
     *                  V(IV,N), W(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, IR, J, JR, LL
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DSPMV
C     .. Executable Statements ..
      IF (NS.EQ.0) RETURN
      LL = 1
      DO 60 I = 1, NS
         IR = IGRADE + I
C
C        FORM THE VECTOR W WHICH CONTAINS THE PRODUCT OF THE
C        IORDER(I) TH COLUMN OF THE MATRIX V AND THE MATRIX H.
C
         CALL DSPMV('Upper triangle',N,1.0D0,H,V(1,IR),1,0.0D0,W,1)
C
C        FORM THE INNER PRODUCT OF W AND THE SEARCH DIRECTION P.
C
         PRHS(I) = DDOT(N,W,1,P,1)
C
C        FORM THE INNER PRODUCT OF W AND THE ROWS OF V.
C
         DO 40 J = 1, I
            JR = IGRADE + J
            SUM = DDOT(N,W,1,V(1,JR),1)
            IF (J.EQ.I) GO TO 20
            PHESL(LL) = SUM
            LL = LL + 1
            GO TO 40
   20       PHESD(I) = SUM
   40    CONTINUE
   60 CONTINUE
      RETURN
C
C     END OF E04HEY   (PRHESS)
C
      END
