      SUBROUTINE C02AFF(A,N,SCALE,Z,WORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15B REVISED. IER-944 (NOV 1991).
C
C     C02AFF ATTEMPTS TO FIND ALL THE ROOTS OF THE NTH ORDER COMPLEX
C     POLYNOMIAL EQUATION
C
C         N
C        SUM [A(1,k)+A(2,k)*I] * Z**(N-k) = 0.
C        k=0
C
C     THE ZEROS OF POLYNOMIALS OF DEGREE 1 AND 2 ARE CALCULATED BY
C     CAREFULLY EVALUATING THE "STANDARD" CLOSED FORMULAS
C           Z = -B/A AND
C           Z = (-B +/- SQRT(B*B-4*A*C))/(2*A) RESPECTIVELY, WHERE
C        A = CMPLX(A(1,0),A(2,0))
C        B = CMPLX(A(1,1),A(2,1)) AND
C        C = CMPLX(A(1,2),A(2,2)).
C     FOR N >= 3, THE ROOTS ARE LOCATED ITERATIVELY USING A VARIANT OF
C     LAGUERRE'S METHOD, WHICH IS CUBICALLY CONVERGENT FOR ISOLATED
C     ZEROS AND LINEARLY CONVERGENT FOR MULTIPLE ZEROS.
C
C     C02AFF ITSELF IS ESSENTIALLY A DUMMY ROUTINE WHOSE FUNCTION IS TO
C     PARTITION THE WORK ARRAY WORK FOR USE BY C02AFZ.
C     WORK IS PARTITIONED INTO 2 ARRAYS EACH OF SIZE 2*(N + 1).
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C02AFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
      LOGICAL           SCALE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(2,0:N), WORK(4*(N+1)), Z(2,N)
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
      EXTERNAL          C02AFZ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IER = IFAIL
      SC = SCALE
      NREC = 0
      NDEG = N
      IF (N.LT.1) THEN
         IER = 1
         WRITE (REC,FMT=99999) N
         NREC = 2
      ELSE IF ((A(1,0).EQ.ZERO) .AND. (A(2,0).EQ.ZERO)) THEN
         IER = 1
         WRITE (REC,FMT=99998)
         NREC = 1
      ELSE
C        Initialize Z to be -infinity.
         BIG = 1.0D0/(SQRT(2.0D0)*X02AMF())
         DO 10 I = 1, N
            Z(1,I) = -BIG
            Z(2,I) = -BIG
   10    CONTINUE
         CALL C02AFZ(A,NDEG,SC,Z,WORK(1),WORK(2*N+3),IER)
         IF (IER.EQ.0) THEN
            IFAIL = 0
            GO TO 20
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
   20 RETURN
C
99999 FORMAT (' ** On entry, N.lt.1:',/'    N = ',I16)
99998 FORMAT (' ** On entry, A(1,0).eq.0 and A(2,0).eq.0')
      END
