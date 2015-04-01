      SUBROUTINE E04UCM(UNITQ,NCQP,NACTIV,NFREE,NZ,N,NLNX,NCTOTL,LDZY,
     *                  LDAQP,LDR,LDT,ISTATE,KACTIV,KX,DXNORM,GDX,AQP,
     *                  ADX,BL,BU,RPQ,RPQ0,DX,GQ,R,T,ZY,WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1082 (JUL 1993).
C     MARK 17 REVISED. IER-1602 (JUN 1995).
C
C     ******************************************************************
C     E04UCM   defines a point which lies on the initial working set for
C     the QP subproblem.  This routine is similar to E04NCH except
C     that advantage is taken of the fact that the initial estimate of
C     the solution of the least-squares subproblem is zero.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level 2 BLAS added 12-June-1986.
C     This version of E04UCM dated 11-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DXNORM, GDX
      INTEGER           LDAQP, LDR, LDT, LDZY, N, NACTIV, NCQP, NCTOTL,
     *                  NFREE, NLNX, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  ADX(*), AQP(LDAQP,*), BL(NCTOTL), BU(NCTOTL),
     *                  DX(N), GQ(N), R(LDR,*), RPQ(NLNX), RPQ0(NLNX),
     *                  T(LDT,*), WORK(N), ZY(LDZY,*)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
      INTEGER           I, J, K, NFIXED, NR
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, DTRMV, E04NBT, E04NBW,
     *                  F06FBF
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
C
      GDX = ZERO
      CALL F06FBF(N,ZERO,DX,1)
      CALL F06FBF(NLNX,ZERO,RPQ,1)
      CALL F06FBF(NLNX,ZERO,RPQ0,1)
C
      IF (NACTIV+NFIXED.GT.0) THEN
C
C        Set  work = residuals for constraints in the working set.
C        Solve for  dx,  the smallest correction to  x  that gives a
C        point on the constraints in the working set.
C        Set the fixed variables on their bounds,  solve the triangular
C        system  T*(dxy) = residuals,  and define  dx = Y*(dxy).
C        Use  (dxy)  to update  d(=Pr)  as  d = d - R'(  0  ).
C                                                     ( dxy )
C
         DO 20 I = 1, NFIXED
            J = KX(NFREE+I)
            IF (ISTATE(J).LE.3) THEN
               BND = BL(J)
               IF (ISTATE(J).EQ.2) BND = BU(J)
               DX(J) = BND
               WORK(NFREE+I) = BND
            ELSE
               WORK(NFREE+I) = ZERO
            END IF
   20    CONTINUE
C
         DO 40 I = 1, NACTIV
            K = KACTIV(I)
            J = N + K
            BND = BL(J)
            IF (ISTATE(J).EQ.2) BND = BU(J)
            WORK(NZ+I) = BND - DDOT(N,AQP(K,1),LDAQP,DX,1)
   40    CONTINUE
C
         IF (NACTIV.GT.0) CALL E04NBT(1,LDT,NACTIV,T(1,NZ+1),WORK(NZ+1))
         CALL DCOPY(NACTIV+NFIXED,WORK(NZ+1),1,DX(NZ+1),1)
         IF (NZ.GT.0) CALL F06FBF(NZ,ZERO,DX,1)
C
         GDX = DDOT(NACTIV+NFIXED,GQ(NZ+1),1,DX(NZ+1),1)
C
         IF (NZ.LT.N) THEN
            CALL DGEMV('N',NZ,N-NZ,-ONE,R(1,NZ+1),LDR,DX(NZ+1),1,ONE,
     *                 RPQ,1)
            IF (NZ.LT.NLNX) THEN
               NR = LDR
               IF (NZ+1.EQ.N) NR = 1
               CALL DCOPY(NLNX-NZ,DX(NZ+1),1,RPQ(NZ+1),1)
               CALL DSCAL(NLNX-NZ,(-ONE),RPQ(NZ+1),1)
               CALL DTRMV('U','N','N',NLNX-NZ,R(NZ+1,NZ+1),NR,RPQ(NZ+1),
     *                    1)
               IF (NLNX.LT.N) THEN
                  NR = LDR
                  IF (NLNX+1.EQ.N) NR = N - NZ
                  CALL DGEMV('N',NLNX-NZ,N-NLNX,-ONE,R(NZ+1,NLNX+1),NR,
     *                       DX(NLNX+1),1,ONE,RPQ(NZ+1),1)
               END IF
            END IF
         END IF
C
         CALL E04NBW(2,N,NZ,NFREE,LDZY,UNITQ,KX,DX,ZY,WORK)
      END IF
C
C     ------------------------------------------------------------------
C     Compute the 2-norm of  DX.
C     Initialize  A*DX.
C     ------------------------------------------------------------------
      DXNORM = DNRM2(N,DX,1)
      IF (NCQP.GT.0) CALL DGEMV('N',NCQP,N,ONE,AQP,LDAQP,DX,1,ZERO,ADX,
     *                          1)
C
      RETURN
C
C
C     End of  E04UCM. (NPSETX)
C
      END
