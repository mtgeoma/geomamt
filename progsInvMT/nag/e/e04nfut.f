      SUBROUTINE E04NFU(N,JTHCOL,H,LDH,X,HX)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04NFU  is used to compute the product Hx, where H is the QP
C     Hessian matrix stored in H and x is an n-vector.
C
C     If LQPTYP is 3 ('QP1') or 4 ('QP2'), the upper-half of the
C     two-dimensional array  H  contains the elements of the
C     symmetric Hessian matrix H.
C     The value of  m  defines the dimension of the m x m leading
C     principal minor of H.
C     The value of m is input as the option 'Hessian Rows' with default
C     value n.  The zero elements of H are not referenced.
C
C     If LQPTYP is 5 ('QP3') or 6 ('QP4'), the Hessian is of the form
C     H = R'R, where  R  is an m x n upper-trapezoidal matrix.  The
C     factor R is stored in the first m rows of the upper half of the
C     two-dimensional array H.  The zero elements of  R  are not
C     accessed.
C
C     Original F66 Version written    March-1982.
C     This version of E04NFU dated 13-May-1990.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           JTHCOL, LDH, N
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LDH,*), HX(N), X(N)
C     .. Scalars in Common ..
      INTEGER           LQPTYP, M
C     .. Local Scalars ..
      INTEGER           JP1, LENRX, NUM
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSYMV, DTRMV, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04NF/LQPTYP, M
C     .. Executable Statements ..
C
      IF (LQPTYP.EQ.3 .OR. LQPTYP.EQ.4) THEN
C
C        Problem type QP1 and QP2.
C
         IF (JTHCOL.GT.0) THEN
C           ------------------------------------------------------------
C           Special case -- extract one column of H.
C           ------------------------------------------------------------
            IF (JTHCOL.GT.M) THEN
               CALL F06FBF(M,(ZERO),HX,1)
            ELSE
               CALL DCOPY(JTHCOL,H(1,JTHCOL),1,HX,1)
               NUM = M - JTHCOL
               JP1 = JTHCOL + 1
               IF (NUM.GT.0) CALL DCOPY(NUM,H(JTHCOL,JP1),LDH,HX(JP1),1)
            END IF
         ELSE
C           ------------------------------------------------------------
C           Normal case.
C           ------------------------------------------------------------
            CALL DSYMV('Upper-triangular',M,ONE,H,LDH,X,1,ZERO,HX,1)
         END IF
C
         IF (N.GT.M) CALL F06FBF(N-M,(ZERO),HX(M+1),1)
C
      ELSE IF (LQPTYP.EQ.5 .OR. LQPTYP.EQ.6) THEN
C
C        Problem type QP3 and QP4.
C
         IF (M.EQ.0) THEN
            CALL F06FBF(N,(ZERO),HX,1)
C
         ELSE IF (JTHCOL.GT.0) THEN
C           ------------------------------------------------------------
C           Special case -- extract one column of H.
C           ------------------------------------------------------------
            LENRX = MIN(JTHCOL,M)
            CALL DCOPY(LENRX,H(1,JTHCOL),1,HX,1)
         ELSE
C           ------------------------------------------------------------
C           Normal case.
C           ------------------------------------------------------------
            CALL DCOPY(M,X,1,HX,1)
            CALL DTRMV('U','N','N',M,H,LDH,HX,1)
            IF (N.GT.M) CALL DGEMV('N',M,N-M,ONE,H(1,M+1),LDH,X(M+1),1,
     *                             ONE,HX,1)
            LENRX = M
         END IF
C
         IF (N.GT.LENRX) CALL DGEMV('T',LENRX,N-LENRX,ONE,H(1,LENRX+1),
     *                              LDH,HX,1,ZERO,HX(LENRX+1),1)
         CALL DTRMV('U','T','N',LENRX,H,LDH,HX,1)
C
      END IF
C
      RETURN
C
C     End of  E04NFU.  (QPHESS)
C
      END
