      SUBROUTINE E04CGZ(N,X,F,G,ISTATE,GPJNRM,COND,POSDEF,NITER,NF,IW,
     *                  LIW,W,LW)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     S-ROUTINE, DUMMY PRINT ROUTINE FOR E04LAF, E04KCF, E04KAF AND
C     E04JAF (I-ROUTINES FOR E04LBR, E04KBN AND E04JBL) AND ALSO FOR
C     E04EBF, E04DFF, E04DEF AND E04CGF.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, F, GPJNRM
      INTEGER           LIW, LW, N, NF, NITER
      LOGICAL           POSDEF
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), W(LW), X(N)
      INTEGER           ISTATE(N), IW(LIW)
C     .. Executable Statements ..
      RETURN
C
C     END OF E04CGZ (SBPRNT)
C
      END
