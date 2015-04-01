      SUBROUTINE G13ASZ(A,LMAX,P,PHI,RR,RI,STAT,INTGR,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LMAX, P
      LOGICAL           STAT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LMAX,LMAX), PHI(P), RI(LMAX), RR(LMAX)
      INTEGER           INTGR(LMAX)
C     .. Local Scalars ..
      INTEGER           I, J, L
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF
      EXTERNAL          A02ABF
C     .. External Subroutines ..
      EXTERNAL          F02AFF
C     .. Executable Statements ..
C
C     initialise A to the zero matrix
C
      DO 40 J = 1, LMAX
         DO 20 I = 1, LMAX
            A(I,J) = 0.0D0
   20    CONTINUE
   40 CONTINUE
C
C     initialise non-zero elements of A
C
      DO 60 L = 1, P
         A(L,1) = PHI(L)
   60 CONTINUE
C
      DO 80 I = 1, LMAX - 1
         A(I,I+1) = 1.0D0
   80 CONTINUE
C
C     calculate eigenvalues of A
C
      IFAIL = 1
      CALL F02AFF(A,LMAX,LMAX,RR,RI,INTGR,IFAIL)
      STAT = .TRUE.
      IF (IFAIL.EQ.1) RETURN
      DO 100 I = 1, LMAX
         IF (A02ABF(RR(I),RI(I)).GE.0.999D0) THEN
            STAT = .FALSE.
            RETURN
         END IF
  100 CONTINUE
      RETURN
C
      END
