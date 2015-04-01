      SUBROUTINE E04MBW(N,NCLIN,NCTOTL,NROWA,LCRASH,LP,MINSUM,NAMED,
     *                  VERTEX,ISTATE,A,AX,BL,BU,CVEC,X)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12B REVISED. IER-537 (FEB 1987).
C
C *********************************************************************
C     E04MBW  PRINTS  A, BL, BU, CVEC, X, A*X,
C                  COLD, LP, MINSUM, NAMED, VERTEX, AND POSSIBLY ISTATE.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF APRIL 1982.  REV. OCT. 1982.
C *********************************************************************
C
C
C     PRINT  WMACH  AND THE LOGICALS.
C
C     .. Scalar Arguments ..
      INTEGER           LCRASH, N, NCLIN, NCTOTL, NROWA
      LOGICAL           LP, MINSUM, NAMED, VERTEX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,N), AX(NROWA), BL(NCTOTL), BU(NCTOTL),
     *                  CVEC(N), X(N)
      INTEGER           ISTATE(NCTOTL)
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  ATX
      INTEGER           I, J, JJ, K, LROWA
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
      CALL X04BAF(NOUT,' ')
      CALL X04BAF(NOUT,' ')
      CALL X04BAF(NOUT,' ')
      CALL X04BAF(NOUT,' ')
      WRITE (REC,FMT=99999)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      CALL X04BAF(NOUT,REC(3))
      DO 20 I = 1, 11
         WRITE (REC,FMT=99998) I, WMACH(I)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
   20 CONTINUE
      WRITE (REC,FMT=99997) LCRASH, LP, MINSUM, NAMED, VERTEX
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
C
C     PRINT  A  BY ROWS AND COMPUTE  AX = A*X.
C
      IF (NCLIN.EQ.0) GO TO 80
      LROWA = NROWA*(N-1) + 1
      DO 60 K = 1, NCLIN
         WRITE (REC,FMT=99996) K
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 40 J = 1, N, 5
            WRITE (REC,FMT=99995) (A(K,JJ),JJ=J,MIN(N,J+4))
            CALL X04BAF(NOUT,REC(1))
   40    CONTINUE
         AX(K) = DDOT(N,A(K,1),NROWA,X,1)
   60 CONTINUE
C
C     PRINT  BL, BU  AND  X OR AX.
C
   80 WRITE (REC,FMT=99994)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 140 J = 1, NCTOTL
         IF (J.GT.N) GO TO 100
         K = J
         ATX = X(J)
         GO TO 120
C
  100    K = J - N
         ATX = AX(K)
         IF (K.EQ.1) THEN
            WRITE (REC,FMT=99993)
            CALL X04BAF(NOUT,REC(1))
            CALL X04BAF(NOUT,REC(2))
         END IF
C
  120    WRITE (REC,FMT=99992) K, BL(J), BU(J), ATX
         CALL X04BAF(NOUT,REC(1))
  140 CONTINUE
C
C     PRINT  CVEC, ISTATE.
C
      IF (LP) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 160 I = 1, N, 5
            WRITE (REC,FMT=99990) (CVEC(J),J=I,MIN(I+4,N))
            CALL X04BAF(NOUT,REC(1))
  160    CONTINUE
      END IF
      IF (LCRASH.GT.0) THEN
         WRITE (REC,FMT=99989)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 180 I = 1, NCTOTL, 10
            WRITE (REC,FMT=99988) (ISTATE(J),J=I,MIN(I+9,NCTOTL))
            CALL X04BAF(NOUT,REC(1))
  180    CONTINUE
      END IF
      RETURN
C
C
C     END OF E04MBW  ( LPDUMP )
99999 FORMAT (' ',/' OUTPUT FROM E04MBW',/' ******************')
99998 FORMAT (/' WMACH(',I2,') =',G15.6)
99997 FORMAT (/' LCRASH =',I3,4X,' LP     =',L3,4X,' MINSUM =',L3,4X,
     *  ' NAMED  =',L3,4X,' VERTEX =',L3)
99996 FORMAT (/' ROW',I6,'  OF  A ...')
99995 FORMAT (5G15.6)
99994 FORMAT (/14X,'J      BL(J)          BU(J)           X(J)')
99993 FORMAT (/14X,'I    BL(N+I)        BU(N+I)         A(I)*X')
99992 FORMAT (I15,3G15.6)
99991 FORMAT (/' CVEC ...')
99990 FORMAT (5G15.6)
99989 FORMAT (/' ISTATE ...')
99988 FORMAT (10I4)
      END
