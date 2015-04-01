      SUBROUTINE G02DCZ(N,ALPHA,X,INCX,A,LDA,C,S,ZETA)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G02DCZ downdates an upper triangular matrix a of order n, for
C     row vector X G02DCZ determineds a orthogonal matrix u  such that
C
C                        (A )     (U  )
C                    Q * (  )  =  (   ) ,
C                        (0 )     ( X )
C
C     where U is upper triangular.
C
C     the matrix Q is determined as the product Q(1)*...*Q(P)
C     where Q(I) is a rotation in the (P+1,I)-plane of the
C     form
C
C                       ( C(I)     -S(I)     )
C                       (                    ) .
C                       ( S(I)       C(I)    )
C
C     The rotations are chosen so that C(I) is real.
C
C
C         If the downdate cannot carried out ZETA is set to
C           -ANORM(INV(A)*X)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, ZETA
      INTEGER           INCX, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), C(N), S(N), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, TEMP
      INTEGER           I, J
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          F06BAF, F06FDF, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
C     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT
C     IN THE ARRAY S.
C
      CALL F06FDF(N,ALPHA,X,INCX,S,1)
      CALL DTRSV('U','T','N',N,A,LDA,S,1)
      ZETA = DNRM2(N,S,1)
      IF (ZETA.LT.1.0D0) THEN
C
C        DETERMINE THE TRANSFORMATIONS.
C
         ZETA = SQRT(1.0D0-ZETA**2)
         TEMP = ZETA
         DO 20 I = N, 1, -1
            AI = S(I)
            CALL F06BAF(TEMP,AI,C(I),S(I))
   20    CONTINUE
C
C        APPLY THE TRANSFORMATIONS TO A.
C
         DO 60 J = 1, N
            AI = 0.0D0
            DO 40 I = J, 1, -1
               TEMP = C(I)*AI + S(I)*A(I,J)
               A(I,J) = C(I)*A(I,J) - S(I)*AI
               AI = TEMP
   40       CONTINUE
   60    CONTINUE
      ELSE
         ZETA = -ZETA
      END IF
      END
