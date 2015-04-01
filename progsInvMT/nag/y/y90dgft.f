      SUBROUTINE Y90DGF(TYPE,NLOSE,A,LDA,M,N)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =============================================================
C         *  Y90DGF :  Chops off binary digits from a complex matrix  *
C         =============================================================
C
C
C     This routines rounds off the least significant bits of the
C     elements of a complex matrix.  The number of bits to be lost is
C     specified in input by the parameter NLOSE.  The whole matrix
C     or each column or row individually can be rounded off (see the
C     input parameter TYPE below).  The largest real/imaginary part of
C     the elements in each column / row / the whole matrix is detected
C     and is rounded off to lose NLOSE bits.  The remaining part of
C     that element and the remaining elementsin the column / row / whole
C     matrix are then rounded to the SAME LEVEL to ensure that all
C     elements in a column / row / the whole matrix have least
C     significant bits of the same order of magnitude in an absolute
C     scale, both for the real and imaginary parts.
C
C
C     Argument List
C     -------------
C
C     NLOSE     Integer, Input.
C               Number of bits to be rounded off with respect to the
C               largest element in each row / column / the whole matrix
C               (according to the value in TYPE).
C
C     A(LDA,*)  Complex*16, Input/Output.
C               The matrix to be rouded off.
C
C     LDA       Integer, Input.
C               Leading dimension of the array A.
C
C     TYPE      Character*1, Input.
C               Defines how the lopping off of the trailing bits is
C               to be carries out.
C               'C'  ==>  Each column individually.
C               'R'  ==>  Each row individually.
C               'A'  ==>  The whole matrix.
C
C     M         Integer, Input.
C               Number of rows of the matrix A.
C
C     N         Integer, Input.
C               Number of columns of the matrix A.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           LDA, M, N, NLOSE
      CHARACTER*1       TYPE
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAX, BASE, FACT, X, Y, Z
      INTEGER           I, J, JEXP, K, L
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, INT, LOG, MAX, MOD, DBLE,
     *                  SIGN
C     .. Statement Functions ..
      DOUBLE PRECISION  FUNC1, FUNC2
C     .. Statement Function definitions ..
      FUNC1(X,Y) = Y + SIGN(X,Y)
      FUNC2(X,Y) = Y + SIGN(X,-Y)
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize
C
C-----------------------------------------------------------------------
      BASE = 2
C-----------------------------------------------------------------------
C
C     Chop each column of the matrix
C
C-----------------------------------------------------------------------
      IF (Y90WAF(TYPE,'C')) THEN
C
         DO 120 J = 1, N
C
C     Calculate the scaling
C
            AMAX = ZERO
            DO 20 I = 1, M
               AMAX = MAX(AMAX,ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
   20       CONTINUE
C
            IF (AMAX.NE.ZERO) THEN
C
               JEXP = INT(LOG(AMAX)/LOG(BASE))
               IF (JEXP.LE.0) JEXP = JEXP - 1
               K = JEXP + NLOSE
               IF (K.GE.0) THEN
                  Z = BASE
               ELSE
                  K = -K
                  Z = ONE/BASE
               END IF
C
               FACT = 1
               DO 40 L = 1, 20
                  IF (K.LE.0) GO TO 60
                  IF (MOD(K,2).GT.0) FACT = FACT*Z
                  Z = Z*Z
                  K = K/2
   40          CONTINUE
   60          CONTINUE
C
C     Now scale the column
C
               IF (NLOSE.GT.0) THEN
                  DO 80 I = 1, M
                     X = FUNC1(FACT,DBLE(A(I,J)))
                     X = FUNC2(FACT,X)
                     Y = FUNC1(FACT,DIMAG(A(I,J)))
                     Y = FUNC2(FACT,Y)
                     A(I,J) = DCMPLX(X,Y)
   80             CONTINUE
               ELSE
                  DO 100 I = 1, M
                     X = DBLE(A(I,J))
                     Y = DIMAG(A(I,J))
                     IF (AMAX/ABS(X).GT.2) THEN
                        X = FUNC1(FACT,X)
                        X = FUNC2(FACT,X)
                     END IF
                     IF (AMAX/ABS(Y).GT.2) THEN
                        Y = FUNC1(FACT,Y)
                        Y = FUNC2(FACT,Y)
                     END IF
                     A(I,J) = DCMPLX(X,Y)
  100             CONTINUE
               END IF
C
            END IF
C
  120    CONTINUE
C-----------------------------------------------------------------------
C
C     Chop each row of the matrix
C
C-----------------------------------------------------------------------
      ELSE IF (Y90WAF(TYPE,'R')) THEN
C
         DO 240 I = 1, M
C
C     Calculate the scaling
C
            AMAX = ZERO
            DO 140 J = 1, N
               AMAX = MAX(AMAX,ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
  140       CONTINUE
C
            IF (AMAX.NE.ZERO) THEN
C
               JEXP = INT(LOG(AMAX)/LOG(BASE))
               IF (JEXP.LE.0) JEXP = JEXP - 1
               K = JEXP + NLOSE
               IF (K.GE.0) THEN
                  Z = BASE
               ELSE
                  K = -K
                  Z = ONE/BASE
               END IF
C
               FACT = 1
               DO 160 L = 1, 20
                  IF (K.LE.0) GO TO 180
                  IF (MOD(K,2).GT.0) FACT = FACT*Z
                  Z = Z*Z
                  K = K/2
  160          CONTINUE
  180          CONTINUE
C
C     Now scale the row
C
               IF (NLOSE.GT.0) THEN
                  DO 200 J = 1, N
                     X = FUNC1(FACT,DBLE(A(I,J)))
                     X = FUNC2(FACT,X)
                     Y = FUNC1(FACT,DIMAG(A(I,J)))
                     Y = FUNC2(FACT,Y)
                     A(I,J) = DCMPLX(X,Y)
  200             CONTINUE
               ELSE
                  DO 220 J = 1, N
                     X = DBLE(A(I,J))
                     Y = DIMAG(A(I,J))
                     IF (AMAX/ABS(X).GT.2) THEN
                        X = FUNC1(FACT,X)
                        X = FUNC2(FACT,X)
                     END IF
                     IF (AMAX/ABS(Y).GT.2) THEN
                        Y = FUNC1(FACT,Y)
                        Y = FUNC2(FACT,Y)
                     END IF
                     A(I,J) = DCMPLX(X,Y)
  220             CONTINUE
               END IF
C
            END IF
C
  240    CONTINUE
C-----------------------------------------------------------------------
C
C     Chop the whole matrix
C
C-----------------------------------------------------------------------
      ELSE
C
C     Calculate the scaling
C
         AMAX = ZERO
         DO 280 J = 1, N
            DO 260 I = 1, M
               AMAX = MAX(AMAX,ABS(DBLE(A(I,J))),ABS(DIMAG(A(I,J))))
  260       CONTINUE
  280    CONTINUE
C
         IF (AMAX.NE.ZERO) THEN
C
            JEXP = INT(LOG(AMAX)/LOG(BASE))
            IF (JEXP.LE.0) JEXP = JEXP - 1
            K = JEXP + NLOSE
            IF (K.GE.0) THEN
               Z = BASE
            ELSE
               K = -K
               Z = ONE/BASE
            END IF
C
            FACT = 1
            DO 300 L = 1, 20
               IF (K.LE.0) GO TO 320
               IF (MOD(K,2).GT.0) FACT = FACT*Z
               Z = Z*Z
               K = K/2
  300       CONTINUE
  320       CONTINUE
C
C     Now scale the matrix
C
            IF (NLOSE.GT.0) THEN
               DO 360 J = 1, N
                  DO 340 I = 1, M
                     X = FUNC1(FACT,DBLE(A(I,J)))
                     X = FUNC2(FACT,X)
                     Y = FUNC1(FACT,DIMAG(A(I,J)))
                     Y = FUNC2(FACT,Y)
                     A(I,J) = DCMPLX(X,Y)
  340             CONTINUE
  360          CONTINUE
            ELSE
               DO 400 J = 1, N
                  DO 380 I = 1, M
                     X = DBLE(A(I,J))
                     Y = DIMAG(A(I,J))
                     IF (AMAX/ABS(X).GT.2) THEN
                        X = FUNC1(FACT,X)
                        X = FUNC2(FACT,X)
                     END IF
                     IF (AMAX/ABS(Y).GT.2) THEN
                        Y = FUNC1(FACT,Y)
                        Y = FUNC2(FACT,Y)
                     END IF
                     A(I,J) = DCMPLX(X,Y)
  380             CONTINUE
  400          CONTINUE
            END IF
C
         END IF
C
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90DGF
C
C-----------------------------------------------------------------------
      RETURN
      END
