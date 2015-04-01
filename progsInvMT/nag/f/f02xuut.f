      SUBROUTINE F02XUU(JOB,D1,E1,DNM1,DN,ENM2,ENM1,C,S)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  F02XUU  determines  the  shift parameters,  c and s,  for  the  first
C  rotation to be performed in a QR step applied to a bidiagonal matrix.
C
C  c and s  are such that
C
C     (  c  s )*( d1**2 - shift ) := ( x )
C     ( -s  c ) (     d1*e1     )    ( 0 )
C
C  and  shift  is  determined  by the  value  of the  parameter  job  as
C  follows:
C
C     When  JOB =  0  then
C
C        shift = 0.
C
C     When  JOB =  1  then
C
C        shift  is the  eigenvalue closest  to the element  a( 2, 2 ) of
C        the matrix A given by
C
C           A = ( dnm1**2 + enm2**2    dnm1*enm1     ).
C               (     dnm1*enm1      dn**2 + enm1**2 )
C
C     When  JOB = -1  then
C
C        shift  is the eigenvalue closest to the element  a( 2, 2 ) of
C        the matrix A given by
C
C           A = ( dnm1**2 + enm1**2  dn*enm1 ).
C               (       dn*enm1       dn**2  )
C
C     When  JOB =  2  then
C
C        shift  is the  smallest eigenvalue  of the  matrix  A  given by
C
C           A = ( dnm1**2 + enm2**2    dnm1*enm1     ).
C               (     dnm1*enm1      dn**2 + enm1**2 )
C
C     When  JOB = -2  then
C
C        shift  is the  smallest eigenvalue  of the  matrix  A  given by
C
C           A = ( dnm1**2 + enm1**2  dn*enm1 ).
C               (       dn*enm1       dn**2  )
C
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 14-September-1987.
C     Sven Hammarling.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, D1, DN, DNM1, E1, ENM1, ENM2, S
      INTEGER           JOB
C     .. Local Scalars ..
      DOUBLE PRECISION  A, APLUSC, B, BOT, CORR, F, RSQEPS, SHIFT,
     *                  SQTEPS, TOP
      LOGICAL           FAIL, FIRST
C     .. External Functions ..
      DOUBLE PRECISION  F06BLF, X02AJF
      EXTERNAL          F06BLF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Save statement ..
      SAVE              FIRST, SQTEPS, RSQEPS
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
      IF (FIRST) THEN
         FIRST = .FALSE.
         SQTEPS = SQRT(X02AJF())
         RSQEPS = 1/SQTEPS
      END IF
C
      IF (JOB.EQ.0) THEN
         A = D1
         B = E1
      ELSE IF (JOB.EQ.1) THEN
         IF (ENM1.EQ.ZERO) THEN
            CORR = ZERO
         ELSE
            TOP = (DNM1-DN)*(DNM1+DN) + (ENM2-ENM1)*(ENM2+ENM1)
            IF (TOP.EQ.ZERO) THEN
               CORR = ENM1*(DNM1-ENM1)
            ELSE
               BOT = 2*DNM1*ENM1
               F = F06BLF(TOP,BOT,FAIL)
               IF (FAIL) THEN
                  CORR = -ENM1**2
               ELSE
                  IF (ABS(F).LT.SQTEPS) THEN
                     BOT = F + SIGN(ONE,F)
                  ELSE IF (ABS(F).GT.RSQEPS) THEN
                     BOT = 2*F
                  ELSE
                     BOT = F + SIGN(ONE,F)*SQRT(1+F**2)
                  END IF
                  CORR = ENM1*(DNM1/BOT-ENM1)
               END IF
            END IF
         END IF
         IF (D1.NE.ZERO) THEN
            A = (ONE-DN/D1)*(D1+DN) + CORR/D1
            B = E1
         ELSE
            A = ONE
            B = ZERO
         END IF
      ELSE IF (JOB.EQ.(-1)) THEN
         TOP = (DN*ENM1)**2
         IF (TOP.EQ.ZERO) THEN
            CORR = ZERO
         ELSE
            F = ((DNM1-DN)*(DNM1+DN)+ENM1**2)/2
            BOT = F + SIGN(ONE,F)*SQRT(TOP+F**2)
            CORR = F06BLF(TOP,BOT,FAIL)
         END IF
         IF (D1.NE.ZERO) THEN
            A = (ONE-DN/D1)*(D1+DN) + CORR/D1
            B = E1
         ELSE
            A = ONE
            B = ZERO
         END IF
      ELSE IF (JOB.EQ.2) THEN
         TOP = 2*((DN*DNM1)**2+(DN*ENM2)**2+(ENM2*ENM1)**2)
         IF (TOP.EQ.ZERO) THEN
            SHIFT = ZERO
         ELSE
            F = ((DNM1-DN)*(DNM1+DN)+(ENM2-ENM1)*(ENM2+ENM1))**2 +
     *          (2*DNM1*ENM1)**2
            APLUSC = DNM1**2 + DN**2 + ENM2**2 + ENM1**2
            BOT = APLUSC + SIGN(ONE,APLUSC)*SQRT(F)
            SHIFT = F06BLF(TOP,BOT,FAIL)
         END IF
         IF (D1.NE.ZERO) THEN
            A = D1 - SHIFT/D1
            B = E1
         ELSE
            A = ONE
            B = ZERO
         END IF
      ELSE IF (JOB.EQ.(-2)) THEN
         TOP = 2*((DN*DNM1)**2)
         IF (TOP.EQ.ZERO) THEN
            SHIFT = ZERO
         ELSE
            F = ((DNM1-DN)**2+ENM1**2)*((DNM1+DN)**2+ENM1**2)
            APLUSC = DNM1**2 + DN**2 + ENM1**2
            BOT = APLUSC + SIGN(ONE,APLUSC)*SQRT(F)
            SHIFT = F06BLF(TOP,BOT,FAIL)
         END IF
         IF (D1.NE.ZERO) THEN
            A = D1 - SHIFT/D1
            B = E1
         ELSE
            A = ONE
            B = ZERO
         END IF
      END IF
      CALL F06BAF(A,B,C,S)
C
      RETURN
C
C     End of F02XUU. ( SBISFT )
C
      END
