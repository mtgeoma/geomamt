      SUBROUTINE G02BXF(WEIGHT,N,M,X,LDX,WT,XBAR,STD,V,LDV,R,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1657 (JUN 1995).
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BXF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDV, LDX, M, N
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  R(LDV,M), STD(M), V(LDV,M), WT(*), X(LDX,M),
     *                  XBAR(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, SW
      INTEGER           I, IERR, IERROR, K, MM, NREC
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FDF, G02AAW, G02BUF, G02BWF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. Executable Statements ..
C
C     INPUT:
C     WEIGHT - Unweighted 'U' , weighted for general weights 'W', for
C              weights 1/var  'V'
C     X      - Data matrix
C     WT     - Vector of weights, optional
C
C     OUTPUT:
C     XBAR   - Mean
C     STD    - Standard deviation
C     V      - Variance matrix
C     R      - Correlation matrix
C
      NREC = 1
      IERROR = 1
      IF (N.LE.1) THEN
         WRITE (REC,FMT=99999) N
      ELSE IF (M.LT.1) THEN
         WRITE (REC,FMT=99998) M
      ELSE IF (LDX.LT.N) THEN
         WRITE (REC,FMT=99997) LDX, N
      ELSE IF (LDV.LT.M) THEN
         WRITE (REC,FMT=99996) LDV, M
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         IERR = 1
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            CALL G02BUF('M','W',N,M,X,LDX,WT,SW,XBAR,R,IERR)
         ELSE IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
            CALL G02BUF('M','U',N,M,X,LDX,WT,SW,XBAR,R,IERR)
         ELSE IF (WEIGHT.EQ.'V' .OR. WEIGHT.EQ.'v') THEN
            CALL G02BUF('M','W',N,M,X,LDX,WT,SW,XBAR,R,IERR)
         ELSE
            IERROR = 2
            WRITE (REC,FMT=99995) WEIGHT
            GO TO 60
         END IF
         IF (IERR.EQ.4) THEN
            IERROR = 3
            WRITE (REC,FMT=99994)
         ELSE
            MM = (M*M+M)/2
            IF (WEIGHT.EQ.'V' .OR. WEIGHT.EQ.'v') THEN
               K = 0
               DO 20 I = 1, N
                  IF (WT(I).GT.0.0D0) K = K + 1
   20          CONTINUE
               IF (K.LE.1) THEN
                  IERROR = 4
                  WRITE (REC,FMT=99991) K
                  GO TO 60
               END IF
               ALPHA = ONE/(DBLE(K)-ONE)
            ELSE
               IF (SW.LE.ONE) THEN
                  IERROR = 4
                  WRITE (REC,FMT=99993) SW
                  GO TO 60
               END IF
               ALPHA = ONE/(SW-ONE)
            END IF
            CALL F06FDF(MM,ALPHA,R,1,V,1)
            CALL G02AAW(M,LDV,V)
            DO 40 I = 1, M
               STD(I) = SQRT(V(I,I))
   40       CONTINUE
            IERR = 1
            CALL G02BWF(M,R,IERR)
            IF (IERR.EQ.2) THEN
               IERROR = 5
               WRITE (REC,FMT=99992)
            END IF
            CALL G02AAW(M,LDV,R)
         END IF
      END IF
C
   60 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
99999 FORMAT (' ** On entry, N.le.1 : N = ',I16)
99998 FORMAT (' ** On entry, M.lt.1 : M = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N : LDX = ',I16,'  N = ',I16)
99996 FORMAT (' ** On entry, LDV.lt.M : LDV = ',I16,'  M = ',I16)
99995 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99994 FORMAT (' ** On entry, there is at least one value of WT.lt.0.0 .'
     *       )
99993 FORMAT (' ** Sum of weights is not greater than 1.0 : SW = ',
     *       D13.5)
99992 FORMAT (' ** On entry, a variable has zero ''variance''.')
99991 FORMAT (' ** Only ',I1,' observation has non-zero weight')
      END
