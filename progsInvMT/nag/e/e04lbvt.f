      SUBROUTINE E04LBV(N,ISTATE,PA,SIGNUM,P)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBV (XPNDPA) EXPANDS THE CORRECTLY ORDERED ELEMENTS OF THE
C     PROJECTED SEARCH DIRECTION INTO THE N-DIMENSIONAL VECTOR
C     NEEDED FOR THE LINEAR SEARCH.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGNUM
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  P(N), PA(N)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      INTEGER           IP, IPA
C     .. Executable Statements ..
      DO 40 IP = 1, N
         IPA = ISTATE(IP)
         IF (IPA.LT.0) GO TO 20
         P(IP) = SIGNUM*PA(IPA)
         GO TO 40
   20    P(IP) = 0.0D+0
   40 CONTINUE
      RETURN
C
C     END OF E04LBV (XPNDPA)
C
      END
