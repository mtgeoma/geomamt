      SUBROUTINE E04NCW(LDH,N,NRANK,TOLRNK,KX,H,MSGLVL,INFORM)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1073 (JUL 1993).
C
C     ******************************************************************
C     E04NCW  forms the Cholesky factorization of the positive
C     semi-definite matrix H such that
C                   PHP'  =  R'R
C     where  P  is a permutation matrix and  R  is upper triangular.
C     The permutation P is chosen to maximize the diagonal of R at each
C     stage.  Only the diagonal and super-diagonal elements of H are
C     used.
C
C     Output:
C
C         INFORM = 0   the factorization was computed successfully,
C                      with the Cholesky factor written in the upper
C                      triangular part of H and P stored in KX.
C                  1   the matrix H was indefinite.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version of E04NCW dated  2-February-1981.
C     Level 2 Blas added 29-June-1986.
C     This version of E04NCW dated  7-July-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOLRNK
      INTEGER           INFORM, LDH, MSGLVL, N, NRANK
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LDH,*)
      INTEGER           KX(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DMAX, SUPMAX
      INTEGER           I, J, K, KMAX
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSWAP, DSYR, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
      INFORM = 0
      NRANK = 0
C
C     Main loop for computing rows of  R.
C
      DO 20 J = 1, N
C
C        Find maximum available diagonal.
C
         KMAX = J - 1 + IDAMAX(N-J+1,H(J,J),LDH+1)
         DMAX = H(KMAX,KMAX)
C
         IF (DMAX.LE.TOLRNK*ABS(H(1,1))) GO TO 40
C
C        Perform a symmetric interchange if necessary.
C
         IF (KMAX.NE.J) THEN
            K = KX(KMAX)
            KX(KMAX) = KX(J)
            KX(J) = K
C
            CALL DSWAP(KMAX-J,H(J+1,KMAX),1,H(J,J+1),LDH)
            CALL DSWAP(J,H(1,J),1,H(1,KMAX),1)
            CALL DSWAP(N-KMAX+1,H(KMAX,KMAX),LDH,H(J,KMAX),LDH)
C
         END IF
C
C        Set the diagonal of  R.
C
         D = SQRT(DMAX)
         H(J,J) = D
         NRANK = NRANK + 1
C
         IF (J.LT.N) THEN
C
C           Set the super-diagonal elements of this row of  R
C           and update the elements of all remaining rows.
C
            CALL DSCAL(N-J,(ONE/D),H(J,J+1),LDH)
            CALL DSYR('U',N-J,-ONE,H(J,J+1),LDH,H(J+1,J+1),LDH)
         END IF
   20 CONTINUE
C
C     Check for the semi-definite case.
C
   40 IF (NRANK.LT.N) THEN
C
C        Find the largest element in the unfactorized block.
C
         SUPMAX = ZERO
         DO 60 I = J, N - 1
            K = I + IDAMAX(N-I,H(I,I+1),LDH)
            SUPMAX = MAX(SUPMAX,ABS(H(I,K)))
   60    CONTINUE
C
         IF (SUPMAX.GT.TOLRNK*ABS(H(1,1))) THEN
            IF (MSGLVL.GT.0) THEN
               WRITE (REC,FMT=99999) DMAX, SUPMAX
               CALL X04BAY(IPRINT,3,REC)
            END IF
            INFORM = 1
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NCW. (LSCHOL)
C
99999 FORMAT (' XXX  Hessian appears to be indefinite.',/' XXX  Maximu',
     *       'm diagonal and off-diagonal ignored in the Cholesky fact',
     *       'orization:',/1X,1P,2D22.14)
      END
