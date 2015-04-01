      SUBROUTINE E04DGT(MSGCG,ITPRT,SLPRT,ITER,NFEVAL,ALFA,OBJF,GNORM,
     *                  OBJGRD,X,XNORM,XKXKNM,N)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-1055 (JUL 1993).
C
C     ******************************************************************
C     This routine prints the various levels of output for E04DGZ.
C     MSGCG controls printout as follows:-
C
C       MSGCG =  0  No printout.
C       MSGCG =  1  Final solution only.
C       MSGCG =  5  One line of output per iteration.
C       MSGCG = 10  Final solution and one line of output per iteration.
C     ******************************************************************
C
C     -- Written on 13th June 1986.
C     Janet Welding, NAG Central Office.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, GNORM, OBJF, XKXKNM, XNORM
      INTEGER           ITER, MSGCG, N, NFEVAL
      LOGICAL           ITPRT, SLPRT
C     .. Array Arguments ..
      DOUBLE PRECISION  OBJGRD(N), X(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      INTEGER           J
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      IF (MSGCG.GE.5 .AND. ITPRT) THEN
         IF (ITER.EQ.0) THEN
C
C           Print heading of iteration print out.
C
            WRITE (REC,FMT=99999)
            CALL X04BAY(NOUT,3,REC)
C
C           STEP and XKXKNM are undefined when ITER = 0.
C
            WRITE (REC,FMT=99998) ITER, NFEVAL, OBJF, GNORM, XNORM
         ELSE
C
C           Print one line per iteration.
C
            WRITE (REC,FMT=99997) ITER, ALFA, NFEVAL, OBJF, GNORM,
     *        XNORM, XKXKNM
         END IF
         CALL X04BAF(NOUT,REC(1))
      ELSE IF (MSGCG.GE.1 .AND. SLPRT) THEN
C
C        Print out heading and final solution.
C
         IF (MSGCG.GE.1) THEN
            WRITE (REC,FMT=99996) ITER
            CALL X04BAY(NOUT,2,REC)
         END IF
         WRITE (REC,FMT=99995)
         CALL X04BAY(NOUT,3,REC)
         DO 20 J = 1, N
            WRITE (REC,FMT=99994) J, X(J), OBJGRD(J)
            CALL X04BAF(NOUT,REC(1))
   20    CONTINUE
      END IF
      SLPRT = .FALSE.
      ITPRT = .FALSE.
      RETURN
C
C     End of E04DGT (CGPRT).
C
99999 FORMAT (//'  Itn      Step  Nfun      Objective    Norm G    Nor',
     *       'm X   Norm (X(k-1)-X(k))')
99998 FORMAT (1X,I4,11X,I5,1X,1P,D14.6,1X,1P,D9.1,1X,1P,D9.1)
99997 FORMAT (1X,I4,1X,1P,D9.1,1X,I5,1X,1P,D14.6,1X,1P,D9.1,1X,1P,D9.1,
     *       7X,1P,D9.1)
99996 FORMAT (/' Exit from E04DGF after ',I5,' iterations.')
99995 FORMAT (//' Variable          Value      Gradient value')
99994 FORMAT (' Varbl',I5,3X,G14.6,4X,1P,D9.1)
      END
