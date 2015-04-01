      SUBROUTINE G11SAS(N2,A,C,LP,P,RL,S,G,ITN,GPJNRM,N3)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     MONITOR E-M ALGORITHM
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GPJNRM
      INTEGER           ITN, LP, N2, N3, S
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N2), C(N2), G(N3), P(S)
      INTEGER           RL(S)
C     .. Local Scalars ..
      DOUBLE PRECISION  LOGL
      INTEGER           I, J, L
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, MIN, DBLE
C     .. Executable Statements ..
      LOGL = 0.0D0
      DO 20 L = 1, S
         LOGL = LOGL + DBLE(RL(L))*LOG(P(L))
   20 CONTINUE
      LOGL = -LOGL
C
      WRITE (REC,FMT=99999) ITN
      CALL X04BAF(LP,REC(1))
      WRITE (REC,FMT=99998) LOGL, GPJNRM
      CALL X04BAY(LP,5,REC)
C
      WRITE (REC,FMT=99997)
      CALL X04BAY(LP,4,REC)
      DO 40 I = 1, N2, 5
         WRITE (REC,FMT=99996) (A(J),J=I,MIN(I+4,N2))
         CALL X04BAF(LP,REC(1))
   40 CONTINUE
      WRITE (REC,FMT=99995)
      CALL X04BAY(LP,4,REC)
      DO 60 I = 1, N2, 5
         WRITE (REC,FMT=99994) (G(2*J-1),J=I,MIN(I+4,N2))
         CALL X04BAF(LP,REC(1))
   60 CONTINUE
      WRITE (REC,FMT=99993)
      CALL X04BAY(LP,4,REC)
      DO 80 I = 1, N2, 5
         WRITE (REC,FMT=99996) (C(J),J=I,MIN(I+4,N2))
         CALL X04BAF(LP,REC(1))
   80 CONTINUE
      WRITE (REC,FMT=99995)
      CALL X04BAY(LP,4,REC)
      DO 100 I = 1, N2, 5
         WRITE (REC,FMT=99994) (G(2*J),J=I,MIN(I+4,N2))
         CALL X04BAF(LP,REC(1))
  100 CONTINUE
C
      WRITE (REC,FMT=99992)
      CALL X04BAY(LP,3,REC)
C
      RETURN
C
99999 FORMAT (' ITERATION NUMBER = ',I5)
99998 FORMAT (/' VALUE OF LOG LIKELIHOOD KERNEL = ',D14.5,//' MAGNITUD',
     *  'E OF LARGEST COMPONENT OF GRADIENT VECTOR = ',D14.3,/)
99997 FORMAT (/' CURRENT ESTIMATES OF ALPHA(J,1)''S',/' --------------',
     *  '-------------------',/)
99996 FORMAT (' ',5F14.3)
99995 FORMAT (/' COMPONENTS OF GRADIENT VECTOR',/' -------------------',
     *  '----------',/)
99994 FORMAT (' ',5D14.3)
99993 FORMAT (/' CURRENT ESTIMATES OF ALPHA(J,0)''S',/' --------------',
     *  '-------------------',/)
99992 FORMAT (/' *****************************************************',
     *  '*******************',/)
      END
