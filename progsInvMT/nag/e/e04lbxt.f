      SUBROUTINE E04LBX(N,ISTATE,GFULL,HDFULL,HLFULL,NH,GPROJ,HDPROJ,
     *                  HLPROJ,LH,LHPROJ)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBX (PACK) PACKS THE RELEVANT ELEMENTS OF THE FULL GRADIENT
C     AND HESSIAN VECTORS IN THE PROJECTED GRADIENT AND HESSIAN
C     VECTORS AND FOR EACH FREE VARIABLE X(IFULL) ENTERS THE
C     POSITIVE INDEX IPROJ IN ISTATE(IFULL). LHPROJ RETURNS THE
C     NUMBER OF ENTRIES IN HLPROJ.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           LH, LHPROJ, N, NH
C     .. Array Arguments ..
      DOUBLE PRECISION  GFULL(N), GPROJ(N), HDFULL(N), HDPROJ(N),
     *                  HLFULL(NH), HLPROJ(LH)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      INTEGER           IFULL, IHLFL, IHLPJ, IMINUS, IPROJ, J
C     .. Executable Statements ..
      IPROJ = 0
      IHLFL = 0
      IHLPJ = 0
      DO 60 IFULL = 1, N
         IMINUS = IFULL - 1
         IF (ISTATE(IFULL).LT.0) GO TO 40
         IPROJ = IPROJ + 1
         ISTATE(IFULL) = IPROJ
         GPROJ(IPROJ) = GFULL(IFULL)
         HDPROJ(IPROJ) = HDFULL(IFULL)
         IF (IFULL.EQ.1) GO TO 60
         DO 20 J = 1, IMINUS
            IHLFL = IHLFL + 1
            IF (ISTATE(J).LT.0) GO TO 20
            IHLPJ = IHLPJ + 1
            HLPROJ(IHLPJ) = HLFULL(IHLFL)
   20    CONTINUE
         GO TO 60
   40    IHLFL = IHLFL + IMINUS
   60 CONTINUE
      LHPROJ = IHLPJ
      RETURN
C
C     END OF E04LBX (PACK)
C
      END
