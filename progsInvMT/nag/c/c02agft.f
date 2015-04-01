      SUBROUTINE C02AGF(A,N,SCALE,Z,WORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15B REVISED. IER-945 (NOV 1991).
C
C     C02AGF ATTEMPTS TO FIND ALL THE ROOTS OF THE NTH ORDER REAL
C     POLYNOMIAL EQUATION
C
C        A(0)*Z**N + A(1)*Z**(N-1) + ... + A(N-1)*Z + A(N) = 0.
C
C     THE ZEROS OF POLYNOMIALS OF DEGREE 1 AND 2 ARE CALCULATED USING
C     THE "STANDARD" CLOSED FORMULAS
C           Z = -A(1)/A(0) AND
C           Z = (-A(1) +/- SQRT(DISC))/(2*A(0)) RESPECTIVELY, WHERE
C        DISC = A(1)**2 - 4*A(0)*A(2).
C     FOR N >= 3, THE ROOTS ARE LOCATED ITERATIVELY USING A VARIANT OF
C     LAGUERRE'S METHOD, WHICH IS CUBICALLY CONVERGENT FOR ISOLATED
C     ZEROS (REAL OR COMPLEX) AND LINEARLY CONVERGENT FOR MULTIPLE
C     ZEROS.
C
C     C02AGF ITSELF IS ESSENTIALLY A DUMMY ROUTINE WHOSE FUNCTION IS TO
C     PARTITION THE WORK ARRAY WORK FOR USE BY C02AGZ.
C     WORK IS PARTITIONED INTO 2 ARRAYS EACH OF SIZE (N + 1).
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C02AGF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
      LOGICAL           SCALE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:N), WORK(2*(N+1)), Z(2,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG
      INTEGER           I, IER, NDEG, NREC
      LOGICAL           SC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C02AGZ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IER = IFAIL
      SC = SCALE
      NREC = 0
      NDEG = N
      IF (N.LT.1 .OR. A(0).EQ.ZERO) THEN
         IER = 1
         WRITE (REC,FMT=99999) N, A(0)
         NREC = 2
         GO TO 20
      END IF
C     Initialize Z to be -infinity.
      BIG = 1.0D0/(SQRT(2.0D0)*X02AMF())
      DO 10 I = 1, N
         Z(1,I) = -BIG
         Z(2,I) = -BIG
   10 CONTINUE
      CALL C02AGZ(A,NDEG,SC,Z,WORK(1),WORK(N+2),IER)
      IF (IER.EQ.0) THEN
         IFAIL = 0
         GO TO 40
      END IF
   20 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
   40 RETURN
C
99999 FORMAT (' ** On entry, N.lt.1 or A(0).eq.0:',/'    N = ',I8,'  A',
     *       '(0) = ',1P,D13.5)
      END
