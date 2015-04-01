      SUBROUTINE F08PSW(SIDE,M,N,V,TAU,C,LDC,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARFX(SIDE,M,N,V,TAU,C,LDC,WORK)
C
C  Purpose
C  =======
C
C  ZLARFX applies a complex elementary reflector H to a complex m by n
C  matrix C, from either the left or the right. H is represented in the
C  form
C
C        H = I - tau * v * v'
C
C  where tau is a complex scalar and v is a complex vector.
C
C  If tau = 0, then H is taken to be the unit matrix
C
C  This version uses inline code if H has order < 11.
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': form  H * C
C          = 'R': form  C * H
C
C  M       (input) INTEGER
C          The number of rows of the matrix C.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C.
C
C  V       (input) COMPLEX*16 array, dimension (M) if SIDE = 'L'
C                                        or (N) if SIDE = 'R'
C          The vector v in the representation of H.
C
C  TAU     (input) COMPLEX*16
C          The value tau in the representation of H.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the m by n matrix C.
C          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
C          or C * H if SIDE = 'R'.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDA >= max(1,M).
C
C  WORK    (workspace) COMPLEX*16 array, dimension (N) if SIDE = 'L'
C                                            or (M) if SIDE = 'R'
C          WORK is not referenced if H has order < 11.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      COMPLEX*16        TAU
      INTEGER           LDC, M, N
      CHARACTER         SIDE
C     .. Array Arguments ..
      COMPLEX*16        C(LDC,*), V(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9,
     *                  V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
      INTEGER           J
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERC
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C
      IF (TAU.EQ.ZERO) RETURN
      IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C        Form  H * C, where H has order m.
C
         GO TO (20,60,100,140,180,220,260,300,340,
     *          380) M
C
C        Code for general M
C
C        w := C'*v
C
         CALL ZGEMV('Conjugate transpose',M,N,ONE,C,LDC,V,1,ZERO,WORK,1)
C
C        C := C - tau * v * w'
C
         CALL ZGERC(M,N,-TAU,V,1,WORK,1,C,LDC)
         GO TO 820
   20    CONTINUE
C
C        Special code for 1 x 1 Householder
C
         T1 = ONE - TAU*V(1)*DCONJG(V(1))
         DO 40 J = 1, N
            C(1,J) = T1*C(1,J)
   40    CONTINUE
         GO TO 820
   60    CONTINUE
C
C        Special code for 2 x 2 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         DO 80 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
   80    CONTINUE
         GO TO 820
  100    CONTINUE
C
C        Special code for 3 x 3 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         DO 120 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
  120    CONTINUE
         GO TO 820
  140    CONTINUE
C
C        Special code for 4 x 4 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         V4 = DCONJG(V(4))
         T4 = TAU*DCONJG(V4)
         DO 160 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J) + V4*C(4,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
            C(4,J) = C(4,J) - SUM*T4
  160    CONTINUE
         GO TO 820
  180    CONTINUE
C
C        Special code for 5 x 5 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         V4 = DCONJG(V(4))
         T4 = TAU*DCONJG(V4)
         V5 = DCONJG(V(5))
         T5 = TAU*DCONJG(V5)
         DO 200 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J) + V4*C(4,J) +
     *            V5*C(5,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
            C(4,J) = C(4,J) - SUM*T4
            C(5,J) = C(5,J) - SUM*T5
  200    CONTINUE
         GO TO 820
  220    CONTINUE
C
C        Special code for 6 x 6 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         V4 = DCONJG(V(4))
         T4 = TAU*DCONJG(V4)
         V5 = DCONJG(V(5))
         T5 = TAU*DCONJG(V5)
         V6 = DCONJG(V(6))
         T6 = TAU*DCONJG(V6)
         DO 240 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J) + V4*C(4,J) +
     *            V5*C(5,J) + V6*C(6,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
            C(4,J) = C(4,J) - SUM*T4
            C(5,J) = C(5,J) - SUM*T5
            C(6,J) = C(6,J) - SUM*T6
  240    CONTINUE
         GO TO 820
  260    CONTINUE
C
C        Special code for 7 x 7 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         V4 = DCONJG(V(4))
         T4 = TAU*DCONJG(V4)
         V5 = DCONJG(V(5))
         T5 = TAU*DCONJG(V5)
         V6 = DCONJG(V(6))
         T6 = TAU*DCONJG(V6)
         V7 = DCONJG(V(7))
         T7 = TAU*DCONJG(V7)
         DO 280 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J) + V4*C(4,J) +
     *            V5*C(5,J) + V6*C(6,J) + V7*C(7,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
            C(4,J) = C(4,J) - SUM*T4
            C(5,J) = C(5,J) - SUM*T5
            C(6,J) = C(6,J) - SUM*T6
            C(7,J) = C(7,J) - SUM*T7
  280    CONTINUE
         GO TO 820
  300    CONTINUE
C
C        Special code for 8 x 8 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         V4 = DCONJG(V(4))
         T4 = TAU*DCONJG(V4)
         V5 = DCONJG(V(5))
         T5 = TAU*DCONJG(V5)
         V6 = DCONJG(V(6))
         T6 = TAU*DCONJG(V6)
         V7 = DCONJG(V(7))
         T7 = TAU*DCONJG(V7)
         V8 = DCONJG(V(8))
         T8 = TAU*DCONJG(V8)
         DO 320 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J) + V4*C(4,J) +
     *            V5*C(5,J) + V6*C(6,J) + V7*C(7,J) + V8*C(8,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
            C(4,J) = C(4,J) - SUM*T4
            C(5,J) = C(5,J) - SUM*T5
            C(6,J) = C(6,J) - SUM*T6
            C(7,J) = C(7,J) - SUM*T7
            C(8,J) = C(8,J) - SUM*T8
  320    CONTINUE
         GO TO 820
  340    CONTINUE
C
C        Special code for 9 x 9 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         V4 = DCONJG(V(4))
         T4 = TAU*DCONJG(V4)
         V5 = DCONJG(V(5))
         T5 = TAU*DCONJG(V5)
         V6 = DCONJG(V(6))
         T6 = TAU*DCONJG(V6)
         V7 = DCONJG(V(7))
         T7 = TAU*DCONJG(V7)
         V8 = DCONJG(V(8))
         T8 = TAU*DCONJG(V8)
         V9 = DCONJG(V(9))
         T9 = TAU*DCONJG(V9)
         DO 360 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J) + V4*C(4,J) +
     *            V5*C(5,J) + V6*C(6,J) + V7*C(7,J) + V8*C(8,J) +
     *            V9*C(9,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
            C(4,J) = C(4,J) - SUM*T4
            C(5,J) = C(5,J) - SUM*T5
            C(6,J) = C(6,J) - SUM*T6
            C(7,J) = C(7,J) - SUM*T7
            C(8,J) = C(8,J) - SUM*T8
            C(9,J) = C(9,J) - SUM*T9
  360    CONTINUE
         GO TO 820
  380    CONTINUE
C
C        Special code for 10 x 10 Householder
C
         V1 = DCONJG(V(1))
         T1 = TAU*DCONJG(V1)
         V2 = DCONJG(V(2))
         T2 = TAU*DCONJG(V2)
         V3 = DCONJG(V(3))
         T3 = TAU*DCONJG(V3)
         V4 = DCONJG(V(4))
         T4 = TAU*DCONJG(V4)
         V5 = DCONJG(V(5))
         T5 = TAU*DCONJG(V5)
         V6 = DCONJG(V(6))
         T6 = TAU*DCONJG(V6)
         V7 = DCONJG(V(7))
         T7 = TAU*DCONJG(V7)
         V8 = DCONJG(V(8))
         T8 = TAU*DCONJG(V8)
         V9 = DCONJG(V(9))
         T9 = TAU*DCONJG(V9)
         V10 = DCONJG(V(10))
         T10 = TAU*DCONJG(V10)
         DO 400 J = 1, N
            SUM = V1*C(1,J) + V2*C(2,J) + V3*C(3,J) + V4*C(4,J) +
     *            V5*C(5,J) + V6*C(6,J) + V7*C(7,J) + V8*C(8,J) +
     *            V9*C(9,J) + V10*C(10,J)
            C(1,J) = C(1,J) - SUM*T1
            C(2,J) = C(2,J) - SUM*T2
            C(3,J) = C(3,J) - SUM*T3
            C(4,J) = C(4,J) - SUM*T4
            C(5,J) = C(5,J) - SUM*T5
            C(6,J) = C(6,J) - SUM*T6
            C(7,J) = C(7,J) - SUM*T7
            C(8,J) = C(8,J) - SUM*T8
            C(9,J) = C(9,J) - SUM*T9
            C(10,J) = C(10,J) - SUM*T10
  400    CONTINUE
         GO TO 820
      ELSE
C
C        Form  C * H, where H has order n.
C
         GO TO (420,460,500,540,580,620,660,700,740,
     *          780) N
C
C        Code for general N
C
C        w := C * v
C
         CALL ZGEMV('No transpose',M,N,ONE,C,LDC,V,1,ZERO,WORK,1)
C
C        C := C - tau * w * v'
C
         CALL ZGERC(M,N,-TAU,WORK,1,V,1,C,LDC)
         GO TO 820
  420    CONTINUE
C
C        Special code for 1 x 1 Householder
C
         T1 = ONE - TAU*V(1)*DCONJG(V(1))
         DO 440 J = 1, M
            C(J,1) = T1*C(J,1)
  440    CONTINUE
         GO TO 820
  460    CONTINUE
C
C        Special code for 2 x 2 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         DO 480 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
  480    CONTINUE
         GO TO 820
  500    CONTINUE
C
C        Special code for 3 x 3 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         DO 520 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
  520    CONTINUE
         GO TO 820
  540    CONTINUE
C
C        Special code for 4 x 4 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         V4 = V(4)
         T4 = TAU*DCONJG(V4)
         DO 560 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3) + V4*C(J,4)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
            C(J,4) = C(J,4) - SUM*T4
  560    CONTINUE
         GO TO 820
  580    CONTINUE
C
C        Special code for 5 x 5 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         V4 = V(4)
         T4 = TAU*DCONJG(V4)
         V5 = V(5)
         T5 = TAU*DCONJG(V5)
         DO 600 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3) + V4*C(J,4) +
     *            V5*C(J,5)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
            C(J,4) = C(J,4) - SUM*T4
            C(J,5) = C(J,5) - SUM*T5
  600    CONTINUE
         GO TO 820
  620    CONTINUE
C
C        Special code for 6 x 6 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         V4 = V(4)
         T4 = TAU*DCONJG(V4)
         V5 = V(5)
         T5 = TAU*DCONJG(V5)
         V6 = V(6)
         T6 = TAU*DCONJG(V6)
         DO 640 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3) + V4*C(J,4) +
     *            V5*C(J,5) + V6*C(J,6)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
            C(J,4) = C(J,4) - SUM*T4
            C(J,5) = C(J,5) - SUM*T5
            C(J,6) = C(J,6) - SUM*T6
  640    CONTINUE
         GO TO 820
  660    CONTINUE
C
C        Special code for 7 x 7 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         V4 = V(4)
         T4 = TAU*DCONJG(V4)
         V5 = V(5)
         T5 = TAU*DCONJG(V5)
         V6 = V(6)
         T6 = TAU*DCONJG(V6)
         V7 = V(7)
         T7 = TAU*DCONJG(V7)
         DO 680 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3) + V4*C(J,4) +
     *            V5*C(J,5) + V6*C(J,6) + V7*C(J,7)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
            C(J,4) = C(J,4) - SUM*T4
            C(J,5) = C(J,5) - SUM*T5
            C(J,6) = C(J,6) - SUM*T6
            C(J,7) = C(J,7) - SUM*T7
  680    CONTINUE
         GO TO 820
  700    CONTINUE
C
C        Special code for 8 x 8 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         V4 = V(4)
         T4 = TAU*DCONJG(V4)
         V5 = V(5)
         T5 = TAU*DCONJG(V5)
         V6 = V(6)
         T6 = TAU*DCONJG(V6)
         V7 = V(7)
         T7 = TAU*DCONJG(V7)
         V8 = V(8)
         T8 = TAU*DCONJG(V8)
         DO 720 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3) + V4*C(J,4) +
     *            V5*C(J,5) + V6*C(J,6) + V7*C(J,7) + V8*C(J,8)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
            C(J,4) = C(J,4) - SUM*T4
            C(J,5) = C(J,5) - SUM*T5
            C(J,6) = C(J,6) - SUM*T6
            C(J,7) = C(J,7) - SUM*T7
            C(J,8) = C(J,8) - SUM*T8
  720    CONTINUE
         GO TO 820
  740    CONTINUE
C
C        Special code for 9 x 9 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         V4 = V(4)
         T4 = TAU*DCONJG(V4)
         V5 = V(5)
         T5 = TAU*DCONJG(V5)
         V6 = V(6)
         T6 = TAU*DCONJG(V6)
         V7 = V(7)
         T7 = TAU*DCONJG(V7)
         V8 = V(8)
         T8 = TAU*DCONJG(V8)
         V9 = V(9)
         T9 = TAU*DCONJG(V9)
         DO 760 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3) + V4*C(J,4) +
     *            V5*C(J,5) + V6*C(J,6) + V7*C(J,7) + V8*C(J,8) +
     *            V9*C(J,9)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
            C(J,4) = C(J,4) - SUM*T4
            C(J,5) = C(J,5) - SUM*T5
            C(J,6) = C(J,6) - SUM*T6
            C(J,7) = C(J,7) - SUM*T7
            C(J,8) = C(J,8) - SUM*T8
            C(J,9) = C(J,9) - SUM*T9
  760    CONTINUE
         GO TO 820
  780    CONTINUE
C
C        Special code for 10 x 10 Householder
C
         V1 = V(1)
         T1 = TAU*DCONJG(V1)
         V2 = V(2)
         T2 = TAU*DCONJG(V2)
         V3 = V(3)
         T3 = TAU*DCONJG(V3)
         V4 = V(4)
         T4 = TAU*DCONJG(V4)
         V5 = V(5)
         T5 = TAU*DCONJG(V5)
         V6 = V(6)
         T6 = TAU*DCONJG(V6)
         V7 = V(7)
         T7 = TAU*DCONJG(V7)
         V8 = V(8)
         T8 = TAU*DCONJG(V8)
         V9 = V(9)
         T9 = TAU*DCONJG(V9)
         V10 = V(10)
         T10 = TAU*DCONJG(V10)
         DO 800 J = 1, M
            SUM = V1*C(J,1) + V2*C(J,2) + V3*C(J,3) + V4*C(J,4) +
     *            V5*C(J,5) + V6*C(J,6) + V7*C(J,7) + V8*C(J,8) +
     *            V9*C(J,9) + V10*C(J,10)
            C(J,1) = C(J,1) - SUM*T1
            C(J,2) = C(J,2) - SUM*T2
            C(J,3) = C(J,3) - SUM*T3
            C(J,4) = C(J,4) - SUM*T4
            C(J,5) = C(J,5) - SUM*T5
            C(J,6) = C(J,6) - SUM*T6
            C(J,7) = C(J,7) - SUM*T7
            C(J,8) = C(J,8) - SUM*T8
            C(J,9) = C(J,9) - SUM*T9
            C(J,10) = C(J,10) - SUM*T10
  800    CONTINUE
         GO TO 820
      END IF
  820 CONTINUE
      RETURN
C
C     End of F08PSW (ZLARFX)
C
      END
