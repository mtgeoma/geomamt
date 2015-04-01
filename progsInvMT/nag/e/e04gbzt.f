      SUBROUTINE E04GBZ(M,N,LH,H,ALPHA,P,GPLUS,G,FJAC,LJ,W,LW)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     **************************************************************
C
C     THIS ROUTINE UPDATES THE SECOND TERM OF THE HESSIAN MATRIX
C     REQUIRED BY THE QUASI-NEWTON LEAST SQUARES ALGORITHM. NOTE
C     THAT THE WORKSPACE ARRAY W MUST BE OF LENGTH AT LEAST
C     M + 3*N.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN AND
C     NICHOLAS I. M. GOULD.
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     Modified to call BLAS.
C     Peter Mayes, NAG Central Office, October 1987.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA
      INTEGER           LH, LJ, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), G(N), GPLUS(N), H(LH), P(N), W(LW)
C     .. Local Scalars ..
      DOUBLE PRECISION  PTQ, YTP
      INTEGER           I, LQ, LY
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DSPMV, DSPR
C     .. Executable Statements ..
C
C     ALLOCATE WORKSPACE ADDRESSES FOR USE WITHIN THE ROUTINE.
C
      LY = M + 1
      LQ = LY + N
C
C     FORM THE VECTOR Y = W(LY) = GK+1 - GK.
C
      DO 20 I = 1, N
         W(LY-1+I) = GPLUS(I) - G(I)
   20 CONTINUE
C
C     FORM THE VECTOR W(LQ) = JT * J * P.
C
      CALL DGEMV('No Transpose',M,N,ONE,FJAC,LJ,P,1,ZERO,W,1)
      CALL DGEMV('Transpose',M,N,ONE,FJAC,LJ,W,1,ZERO,W(LQ),1)
C
C     FORM THE VECTOR Q = W(LQ) = JT * J * P + H * P.
C
      CALL DSPMV('Upper Triangle',N,ONE,H,P,1,ONE,W(LQ),1)
C
C     FORM THE SCALARS ALPHA * YT * P AND PT * Q
C
      YTP = 2.0D0*ALPHA*DDOT(N,W(LY),1,P,1)
      PTQ = DDOT(N,W(LQ),1,P,1)
C
C     UPDATE THE MATRIX H USING THE TRANSFORMATION HK+1 =
C     H + Y * YT/YTP - Q * QT/PTQ, WHERE T DENOTES THE TRANSPOSE.
C
      CALL DSPR('Upper Triangle',N,ONE/YTP,W(LY),1,H)
      CALL DSPR('Upper Triangle',N,-ONE/PTQ,W(LQ),1,H)
      RETURN
C
C     END OF E04GBZ   (UPFHES)
C
      END
