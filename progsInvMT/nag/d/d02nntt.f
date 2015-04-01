      SUBROUTINE D02NNT(N,NZ,LICN,LIRN,IA,JA,IWK,IAMY,A,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C-----------------------------------------------------------------------
C     IA    ROW POINTERS
C     IAMY  COMPACT ROW POINTERS
C     JA    COLUMN POINTERS
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LICN, LIRN, N, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NZ)
      INTEGER           IA(NZ), IAMY(N+1), IWK(N), JA(NZ)
C     .. Scalars in Common ..
      INTEGER           I1, IDSTRC
C     .. Local Scalars ..
      INTEGER           I, ICOL, J, K, KEND
C     .. Common blocks ..
      COMMON            /GD02NN/IDSTRC
      COMMON            /ZD02NN/I1
C     .. Save statement ..
      SAVE              /GD02NN/, /ZD02NN/
C     .. Executable Statements ..
      IF (N.GT.0) GO TO 20
      IFAIL = 1 + (1-IDSTRC)*7
      RETURN
C
   20 IF (NZ.GT.0) GO TO 40
      IFAIL = 2 + (1-IDSTRC)*7
      RETURN
C
   40 IF (LIRN.GE.NZ) GO TO 60
      IFAIL = 3 + (1-IDSTRC)*7
      RETURN
C
   60 IF (LICN.GE.NZ) GO TO 80
      IFAIL = 11
      RETURN
C
C CHECK THAT ROW/COL INDICES ARE WITHIN THE RANGE 1 TO N
C
   80 DO 100 I = 1, NZ
         IF (IA(I).GT.0 .AND. IA(I).LE.N .AND. JA(I).GT.0 .AND. JA(I)
     *       .LE.N) GO TO 100
         I1 = I
         IFAIL = 4 + (1-IDSTRC)*8
         RETURN
  100 CONTINUE
C
C CHECK FOR DUPLICATE ELEMENTS
C
      K = IAMY(1)
      DO 160 I = 1, N
         KEND = IAMY(I+1)
         DO 120 J = 1, N
            IWK(J) = 0
  120    CONTINUE
         DO 140 J = K, KEND - 1
            ICOL = JA(J)
            IF (IWK(ICOL).GT.0) THEN
               I1 = I
               IFAIL = 8 + (1-IDSTRC)*5
               RETURN
            ELSE
               IWK(ICOL) = 1
            END IF
  140    CONTINUE
         K = KEND
  160 CONTINUE
      RETURN
      END
