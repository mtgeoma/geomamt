      SUBROUTINE C05NCZ(M,N,S,LS,U,V,W,SING)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05NCZ (based on MINPACK Routine R1UPDT)
C
C     Given an M by N lower trapezoidal matrix S, an M-vector U,
C     and an N-vector V, the problem is to determine an
C     orthogonal matrix Q such that
C
C             T
C     (S + U*V )*Q
C
C     is again lower trapezoidal.
C
C     This subroutine determines Q as the product of 2*(N - 1)
C     transformations
C
C        GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     where GV(I), GW(I) are Givens rotations in the (I,N) plane
C     which eliminate elements in the I-th and N-th planes,
C     respectively. Q itself is not accumulated, rather the
C     information to recover the GV, GW rotations is returned.
C
C     The subroutine statement is
C
C        SUBROUTINE C05NCZ(M,N,S,LS,U,V,W,SING)
C
C     where
C
C     M is a positive integer input variable set to the number
C     of rows of S.
C
C     N is a positive integer input variable set to the number
C     of columns of S. N must not exceed M.
C
C     S is an array of length LS. On input S must contain the lower
C     trapezoidal matrix S stored by columns. On output S contains
C     the lower trapezoidal matrix produced as described above.
C
C     LS is a positive integer input variable not less than
C     (N*(2*M-N+1))/2.
C
C     U is an input array of length M which must contain the
C     vector U.
C
C     V is an array of length N. On input V must contain the vector
C     V. On output V(I) contains the information necessary to
C     recover the Givens rotation GV(I) described above.
C
C     W is an output array of length M. W(I) contains information
C     necessary to recover the Givens rotation GW(I) described
C     above.
C
C     SING is a logical output variable. SING is set .TRUE. if any
C     of the diagonal elements of the output S are zero. Otherwise
C     SING is set .FALSE.
C
C     Argonne National Laboratory. MINPACK project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More,
C     John L. Nazareth
C
C     **********
C
C     Revised to call BLAS.
C     P.J.D. Mayes, J.J. Du Croz, NAG Central Office, September 1987.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, P5, P25, ZERO
      PARAMETER         (ONE=1.0D0,P5=0.5D0,P25=0.25D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           LS, M, N
      LOGICAL           SING
C     .. Array Arguments ..
      DOUBLE PRECISION  S(LS), U(M), V(N), W(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  COSINE, COTAN, GIANT, SINE, TANGNT, TAU
      INTEGER           I, J, JJ, L, NM1, NMJ
C     .. External Functions ..
      DOUBLE PRECISION  X02ALF
      EXTERNAL          X02ALF
C     .. External Subroutines ..
      EXTERNAL          DROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      GIANT = X02ALF()
C
C     Initialize the diagonal element pointer.
C
      JJ = (N*(2*M-N+1))/2 - (M-N)
C
C     Move the nontrivial part of the last column of S into W.
C
      L = JJ
      DO 20 I = N, M
         W(I) = S(L)
         L = L + 1
   20 CONTINUE
C
C     Rotate the vector V into a multiple of the N-th unit vector
C     in such a way that a spike is introduced into W.
C
      NM1 = N - 1
      IF (NM1.GE.1) THEN
         DO 60 NMJ = 1, NM1
            J = N - NMJ
            JJ = JJ - (M-J+1)
            W(J) = ZERO
            IF (V(J).EQ.ZERO) GO TO 40
C
C           Determine a Givens rotation which eliminates the
C           J-th element of V.
C
            IF (ABS(V(N)).LT.ABS(V(J))) THEN
               COTAN = V(N)/V(J)
               SINE = P5/SQRT(P25+P25*COTAN**2)
               COSINE = SINE*COTAN
               TAU = ONE
               IF (ABS(COSINE)*GIANT.GT.ONE) TAU = ONE/COSINE
            ELSE
               TANGNT = V(J)/V(N)
               COSINE = P5/SQRT(P25+P25*TANGNT**2)
               SINE = COSINE*TANGNT
               TAU = SINE
            END IF
C
C           Apply the transformation to V and store the information
C           necessary to recover the Givens rotation.
C
            V(N) = SINE*V(J) + COSINE*V(N)
            V(J) = TAU
C
C           Apply the transformation to S and extend the spike in W.
C
            CALL DROT(M-J+1,W(J),1,S(JJ),1,COSINE,SINE)
   40       CONTINUE
   60    CONTINUE
      END IF
C
C     Add the spike from the rank 1 update to W.
C
      DO 80 I = 1, M
         W(I) = W(I) + V(N)*U(I)
   80 CONTINUE
C
C     Eliminate the spike.
C
      SING = .FALSE.
      IF (NM1.GE.1) THEN
         DO 100 J = 1, NM1
            IF (W(J).NE.ZERO) THEN
C
C              Determine a Givens rotation which eliminates the
C              J-th element of the spike.
C
               IF (ABS(S(JJ)).LT.ABS(W(J))) THEN
                  COTAN = S(JJ)/W(J)
                  SINE = P5/SQRT(P25+P25*COTAN**2)
                  COSINE = SINE*COTAN
                  TAU = ONE
                  IF (ABS(COSINE)*GIANT.GT.ONE) TAU = ONE/COSINE
               ELSE
                  TANGNT = W(J)/S(JJ)
                  COSINE = P5/SQRT(P25+P25*TANGNT**2)
                  SINE = COSINE*TANGNT
                  TAU = SINE
               END IF
C
C              Apply the transformation to S and reduce the spike in W.
C
               CALL DROT(M-J+1,S(JJ),1,W(J),1,COSINE,SINE)
C
C              Store the information necessary to recover the
C              Givens rotation.
C
               W(J) = TAU
            END IF
C
C           Test for zero diagonal elements in the output S.
C
            IF (S(JJ).EQ.ZERO) SING = .TRUE.
            JJ = JJ + (M-J+1)
  100    CONTINUE
      END IF
C
C     Move W back into the last column of the output S.
C
      L = JJ
      DO 120 I = N, M
         S(L) = W(I)
         L = L + 1
  120 CONTINUE
      IF (S(JJ).EQ.ZERO) SING = .TRUE.
      RETURN
      END
