      SUBROUTINE S17DGS(ZR,S1,S2,NZ,ASCLE,ALIM,IUF)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-772 (DEC 1989).
C
C     Original name: CS1S2
C
C     S17DGS TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C     .. Scalar Arguments ..
      COMPLEX*16        S1, S2, ZR
      DOUBLE PRECISION  ALIM, ASCLE
      INTEGER           IUF, NZ
C     .. Local Scalars ..
      COMPLEX*16        C1, CZERO, S1D
      DOUBLE PRECISION  AA, ALN, AS1, AS2, XX
      INTEGER           IF1
C     .. External Functions ..
      COMPLEX*16        S01EAF
      EXTERNAL          S01EAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, LOG, MAX, DBLE
C     .. Data statements ..
      DATA              CZERO/(0.0D0,0.0D0)/
C     .. Executable Statements ..
C
      NZ = 0
      AS1 = ABS(S1)
      AS2 = ABS(S2)
      AA = DBLE(S1)
      ALN = DIMAG(S1)
      IF (AA.NE.0.0D0 .OR. ALN.NE.0.0D0) THEN
         IF (AS1.NE.0.0D0) THEN
            XX = DBLE(ZR)
            ALN = -XX - XX + LOG(AS1)
            S1D = S1
            S1 = CZERO
            AS1 = 0.0D0
            IF (ALN.GE.(-ALIM)) THEN
               C1 = LOG(S1D) - ZR - ZR
C               S1 = EXP(C1)
               IF1 = 1
               S1 = S01EAF(C1,IF1)
               AS1 = ABS(S1)
               IUF = IUF + 1
            END IF
         END IF
      END IF
      AA = MAX(AS1,AS2)
      IF (AA.LE.ASCLE) THEN
         S1 = CZERO
         S2 = CZERO
         NZ = 1
         IUF = 0
      END IF
      RETURN
      END
