      SUBROUTINE E04JBX(LHD,HESD,N,HESL,NH,COND)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04JBX (UNITHS) SETS THE CHOLESKY FACTORS FOR THE INITIAL
C     APPROXIMATION OF THE PROJECTED HESSIAN TO THOSE OF THE UNIT
C     MATRIX, WITH CONDITION NUMBER 1, AND THE REMAINING ELEMENTS OF
C     HESD AND HESL TO 0.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND
      INTEGER           LHD, N, NH
C     .. Array Arguments ..
      DOUBLE PRECISION  HESD(N), HESL(NH)
C     .. Local Scalars ..
      INTEGER           I, LHD1, LHL
C     .. Executable Statements ..
      IF (LHD.EQ.0) GO TO 40
      DO 20 I = 1, LHD
         HESD(I) = 1.0D+0
   20 CONTINUE
      IF (LHD.EQ.N) GO TO 80
   40 LHD1 = LHD + 1
      DO 60 I = LHD1, N
         HESD(I) = 0.0D+0
   60 CONTINUE
   80 LHL = 1
      IF (N.GE.3) LHL = N*(N-1)/2
      DO 100 I = 1, LHL
         HESL(I) = 0.0D+0
  100 CONTINUE
      COND = 1.0D+0
      RETURN
C
C     END OF E04JBX (UNITHS)
C
      END
