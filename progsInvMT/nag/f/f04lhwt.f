      SUBROUTINE F04LHW(UPLO,TRANS,DIAG,N,M1,M2,A,LDA,X,INCX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Purpose
C     =======
C
C     F04LHW  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C     where b and x are n element vectors and A is an n by n unit, or
C     non-unit, upper or lower triangular matrix of the special forms
C
C     (  I   A12  A13 )  (upper)  or  (  I            )  (lower)
C     (      A22  A23 )               ( A21  A22      )
C     (            I  )               ( A31  A32   I  )
C
C     where the submatrices A22 in rows and columns m1 to m2 are
C     upper or lower triangular and the off-diagonal submatrices
C     are rectangular.
C
C     No test for singularity or near-singularity is included in this
C     routine. Such tests must be performed before calling this routine.
C
C     Parameters
C     ==========
C
C     UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C     TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C     DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C     N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C     M1     - INTEGER
C           On entry, M1 specifies the first row and column of the
C           diagonal submatrix A22. M1 must be at least 1.
C           Unchanged on exit.
C
C     M2     - INTEGER
C           On entry, M2 specifies the last row and column of the
C           diagonal submatrix A22. M2 must be at least m1 - 1,
C           and must not exceed n.
C           Unchanged on exit.
C
C     A      - REAL             array of DIMENSION ( LDA, n ).
C           Before entry with UPLO = 'U' or 'u', rows 1 to m2 and column
C           m1 to n of the upper triangular part of the array A must
C           contain the elements of the submatrices A12, A13, A22 and A2
C           The rest of the array is not referenced.
C           Before entry with UPLO = 'L' or 'l', rows m1 to n and column
C           1 to m2 of the lower triangular part of the array A must
C           contain the elements of the submatrices A21, A22, A31 and A3
C           The rest of the array is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A22 are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C     LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least n.
C           Unchanged on exit.
C
C     X      - REAL             array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C     INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X.
C           Unchanged on exit.
C
C
C     Note that UPLO, TRANS, DIAG, N, M1, M2 and LDA must be such that
C     the value of the LOGICAL variable OK in the following statement
C     is true.
C
C      OK = ( ( UPLO.EQ.'U' ).OR.( UPLO.EQ.'u' ).OR.
C     $       ( UPLO.EQ.'L' ).OR.( UPLO.EQ.'l' )     )
C     $     .AND.
C     $     ( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ).OR.
C     $       ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
C     $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )
C     $     .AND.
C     $     ( ( DIAG.EQ.'U' ).OR.( DIAG.EQ.'u' ).OR.
C     $       ( DIAG.EQ.'N' ).OR.( DIAG.EQ.'n' )     )
C     $     .AND.
C     $     ( N.GE.0 )
C     $     .AND.
C     $     ( LDA.GE.N )
C     $     .AND.
C     $     ( M1.GE.1 )
C     $     .AND.
C     $     ( M2.GE.M1-1 )
C     $     .AND.
C     $     ( M2.LE.N )
C
C     Level 2 Blas routine.
C
C     -- Written on 19-March-1986.
C     Sven Hammarling and Jeremy Du Croz, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INCX, LDA, M1, M2, N
      CHARACTER*1       DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), X(*)
C     .. Local Scalars ..
      INTEGER           I, IX, J, JX, KX
      LOGICAL           NOUNIT
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF (N.EQ.0) RETURN
      NOUNIT = (DIAG.EQ.'N') .OR. (DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF (INCX.LE.0) THEN
         KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF ((TRANS.EQ.'N') .OR. (TRANS.EQ.'n')) THEN
C
C        Form  x := inv( A )*x.
C
         IF ((UPLO.EQ.'U') .OR. (UPLO.EQ.'u')) THEN
            IF (INCX.EQ.1) THEN
               DO 40 J = N, M1, -1
                  IF (X(J).NE.ZERO) THEN
                     IF (NOUNIT .AND. J.LE.M2) X(J) = X(J)/A(J,J)
                     DO 20 I = MIN(M2,J-1), 1, -1
                        X(I) = X(I) - X(J)*A(I,J)
   20                CONTINUE
                  END IF
   40          CONTINUE
            ELSE
               JX = KX + (N-1)*INCX
               KX = KX + M2*INCX
               DO 80 J = N, M1, -1
                  IF (X(JX).NE.ZERO) THEN
                     IF (NOUNIT .AND. J.LE.M2) X(JX) = X(JX)/A(J,J)
                     IX = MIN(KX,JX)
                     DO 60 I = MIN(M2,J-1), 1, -1
                        IX = IX - INCX
                        X(IX) = X(IX) - X(JX)*A(I,J)
   60                CONTINUE
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
               DO 120 J = 1, M2
                  IF (X(J).NE.ZERO) THEN
                     IF (NOUNIT .AND. J.GE.M1) X(J) = X(J)/A(J,J)
                     DO 100 I = MAX(M1,J+1), N
                        X(I) = X(I) - X(J)*A(I,J)
  100                CONTINUE
                  END IF
  120          CONTINUE
            ELSE
               JX = KX
               KX = KX + (M1-2)*INCX
               DO 160 J = 1, M2
                  IF (X(JX).NE.ZERO) THEN
                     IF (NOUNIT .AND. J.GE.M1) X(JX) = X(JX)/A(J,J)
                     IX = MAX(KX,JX)
                     DO 140 I = MAX(M1,J+1), N
                        IX = IX + INCX
                        X(IX) = X(IX) - X(JX)*A(I,J)
  140                CONTINUE
                  END IF
                  JX = JX + INCX
  160          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x.
C
         IF ((UPLO.EQ.'U') .OR. (UPLO.EQ.'u')) THEN
            IF (INCX.EQ.1) THEN
               DO 200 J = M1, N
                  DO 180 I = 1, MIN(M2,J-1)
                     X(J) = X(J) - A(I,J)*X(I)
  180             CONTINUE
                  IF (NOUNIT .AND. J.LE.M2) X(J) = X(J)/A(J,J)
  200          CONTINUE
            ELSE
               JX = KX + (M1-1)*INCX
               DO 240 J = M1, N
                  IX = KX
                  DO 220 I = 1, MIN(M2,J-1)
                     X(JX) = X(JX) - A(I,J)*X(IX)
                     IX = IX + INCX
  220             CONTINUE
                  IF (NOUNIT .AND. J.LE.M2) X(JX) = X(JX)/A(J,J)
                  JX = JX + INCX
  240          CONTINUE
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
               DO 280 J = M2, 1, -1
                  DO 260 I = N, MAX(M1,J+1), -1
                     X(J) = X(J) - A(I,J)*X(I)
  260             CONTINUE
                  IF (NOUNIT .AND. J.GE.M1) X(J) = X(J)/A(J,J)
  280          CONTINUE
            ELSE
               JX = KX + (M2-1)*INCX
               KX = KX + (N-1)*INCX
               DO 320 J = M2, 1, -1
                  IX = KX
                  DO 300 I = N, MAX(M1,J+1), -1
                     X(JX) = X(JX) - A(I,J)*X(IX)
                     IX = IX - INCX
  300             CONTINUE
                  IF (NOUNIT .AND. J.GE.M1) X(JX) = X(JX)/A(J,J)
                  JX = JX - INCX
  320          CONTINUE
            END IF
         END IF
      END IF
      RETURN
C
C     End of F04LHW .
C
      END
