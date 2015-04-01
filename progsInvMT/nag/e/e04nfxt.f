      SUBROUTINE E04NFX(UNITQ,QPHESS,MAXNZ,N,NGQ,NRZ,NZ,NFREE,LDQ,LDH,
     *                  LDR,KX,HSIZE,TOLRNK,GQ,H,R,Q,HZ,WRK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1596 (JUN 1995).
C
C     ******************************************************************
C     E04NFX  computes the Cholesky factor Rz of the reduced Hessian
C     Z'HZ,  given the (NFREE x NZ) matrix Z.  If the reduced Hessian is
C     indefinite, the Cholesky factor of the (NRZ x NRZ) matrix  Hz  is
C     returned, where Hz is formed from  H  and  NRZ  columns of Z.
C     Column interchanges are used in an attempt to maximize NRZ.
C     These are applied to  Z  and the rows of the matrix  GQ.
C
C     This version of E04NFX dated 13-Dec-1990.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HSIZE, TOLRNK
      INTEGER           LDH, LDQ, LDR, MAXNZ, N, NFREE, NGQ, NRZ, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  GQ(N,*), H(LDH,*), HZ(N), Q(LDQ,*), R(LDR,*),
     *                  WRK(N)
      INTEGER           KX(N)
C     .. Subroutine Arguments ..
      EXTERNAL          QPHESS
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DMAX, DMIN
      INTEGER           I, J, JTHCOL, K, KMAX, MNZ
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, DSWAP, DSYR, E04NBW, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      NRZ = 0
      IF (NZ.EQ.0) GO TO 80
C
      MNZ = MIN(NZ,MAXNZ)
C
C     ------------------------------------------------------------------
C     Compute  Z'HZ  and store the upper-triangular symmetric part in
C     the first  MNZ  columns of R.
C     ------------------------------------------------------------------
      DO 40 K = 1, MNZ
         CALL F06FBF(N,(ZERO),WRK,1)
         IF (UNITQ) THEN
C
C           Only bounds are in the working set.  The k-th column of Z
C           is just a column of the identity matrix.
C
            JTHCOL = KX(K)
            WRK(JTHCOL) = ONE
         ELSE
C
C           Expand the column of Z into an n-vector.
C
            DO 20 I = 1, NFREE
               J = KX(I)
               WRK(J) = Q(I,K)
   20       CONTINUE
            JTHCOL = 0
         END IF
C
C        Set  R(*,k)  =  top of  H*(column of Z).
C
         CALL QPHESS(N,JTHCOL,H,LDH,WRK,HZ)
         CALL E04NBW(4,N,NZ,NFREE,LDQ,UNITQ,KX,HZ,Q,WRK)
         CALL DCOPY(MNZ,HZ,1,R(1,K),1)
C
C        Update an estimate of the size of the reduced Hessian.
C
         HSIZE = MAX(HSIZE,ABS(R(K,K)))
   40 CONTINUE
C
C     ------------------------------------------------------------------
C     Form the Cholesky factorization R'R = Z'HZ as far as possible,
C     using symmetric interchanges.
C     ------------------------------------------------------------------
      DMIN = TOLRNK*HSIZE
C
      DO 60 J = 1, MNZ
C
C        Find the maximum diagonal of the Schur complement.
C
         KMAX = J - 1 + IDAMAX(MNZ-J+1,R(J,J),LDR+1)
         DMAX = R(KMAX,KMAX)
C
C        See if the diagonal is big enough.
C
         IF (DMAX.LE.DMIN) GO TO 80
C
C        Perform a symmetric interchange if necessary.
C
         IF (KMAX.NE.J) THEN
            IF (UNITQ) THEN
               K = KX(KMAX)
               KX(KMAX) = KX(J)
               KX(J) = K
            ELSE
               CALL DSWAP(NFREE,Q(1,KMAX),1,Q(1,J),1)
            END IF
C
            IF (NGQ.GT.0) CALL DSWAP(NGQ,GQ(KMAX,1),N,GQ(J,1),N)
            CALL DSWAP(KMAX-J,R(J+1,KMAX),1,R(J,J+1),LDR)
            CALL DSWAP(J,R(1,J),1,R(1,KMAX),1)
            CALL DSWAP(MNZ-KMAX+1,R(KMAX,KMAX),LDR,R(J,KMAX),LDR)
         END IF
C
C        Set the diagonal element of R.
C
         D = SQRT(DMAX)
         R(J,J) = D
         NRZ = NRZ + 1
C
         IF (J.LT.MNZ) THEN
C
C           Set the super-diagonal elements of this row of R and update
C           the elements of the Schur complement.
C
            CALL DSCAL(MNZ-J,(ONE/D),R(J,J+1),LDR)
            CALL DSYR('U',MNZ-J,(-ONE),R(J,J+1),LDR,R(J+1,J+1),LDR)
         END IF
   60 CONTINUE
C
   80 RETURN
C
C
C     End of E04NFX.  (QPCRSH)
C
      END
