      SUBROUTINE Y90RCF(TYPE,DTYPE,M,N,A,IA,D,DIAG,COND,SCALE,DETMAN,
     *                  DETEXP,DIST,SEED,PVROW,PVCOL,IPVROW,IPVCOL,U,IU,
     *                  V,IV,WORK1,IWORK1,WORK2,IWORK2)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =======================================
C         *  Y90RCF :  Random Matrix Generator  *
C         =======================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DETMAN, SCALE
      INTEGER           DETEXP, DIST, DTYPE, IA, IU, IV, IWORK1, IWORK2,
     *                  M, N, TYPE
      CHARACTER*1       PVCOL, PVROW
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), D(*), DIAG(*), U(IU,*), V(IV,*),
     *                  WORK1(IWORK1,*), WORK2(IWORK2,*)
      INTEGER           IPVCOL(*), IPVROW(*), SEED(4)
C     .. Local Scalars ..
      INTEGER           I, J, ND
      CHARACTER*1       DETREQ
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F06EDF, F06EGF, F06YAF, Y90RDF, Y90RGF, Y90SHF,
     *                  Y90SKF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize
C
C-----------------------------------------------------------------------
      ND = MIN(M,N)
C-----------------------------------------------------------------------
C
C     Generate diagonal
C
C-----------------------------------------------------------------------
      IF (TYPE.NE.1 .AND. DTYPE.NE.0) THEN
         IF (M.EQ.N) THEN
            DETREQ = 'D'
         ELSE
            DETREQ = 'N'
         END IF
         CALL Y90RGF(DETREQ,DTYPE,N,D,DIAG,COND,SCALE,DETMAN,DETEXP,
     *               DIST,SEED)
      END IF
C-----------------------------------------------------------------------
C
C     Generate Random or Orthogonal Matrix
C
C-----------------------------------------------------------------------
      IF (TYPE.NE.0) THEN
C
         IF (M.GT.1) THEN
            DO 40 I = 1, M
               DO 20 J = 1, N
                  A(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
         END IF
C
         IF (TYPE.EQ.1) THEN
            IF (M.GE.N) THEN
               CALL Y90RDF('Left','Initialize',M,N,A,IA,WORK1,WORK2,
     *                     SEED)
            ELSE
               CALL Y90RDF('Right','Initialize',M,N,A,IA,WORK1,WORK2,
     *                     SEED)
            END IF
C
         ELSE IF (TYPE.EQ.2) THEN
            DO 60 I = 1, ND
               A(I,I) = D(I)
   60       CONTINUE
            CALL Y90RDF('Left','No Initialize',M,N,A,IA,WORK1,WORK2,
     *                  SEED)
            CALL Y90RDF('Right','No Initialize',M,N,A,IA,WORK1,WORK2,
     *                  SEED)
C
         ELSE IF (TYPE.EQ.3) THEN
            CALL Y90RDF('Left','Initialize',M,N,U,IU,WORK1,WORK2,SEED)
            CALL Y90SHF('T',M,N,U,IU,V,IV)
            DO 80 I = 1, M
               CALL F06EDF(M,D(I),U(1,I),1)
   80       CONTINUE
            CALL F06YAF('N','N',M,M,M,ONE,U,IU,V,IV,ZERO,A,IA)
C
C     Symmetrize the matrix A
C
            CALL Y90SKF('U',A,IA,N)
         END IF
C-----------------------------------------------------------------------
C
C     Permute Rows and/or columns when required
C
C-----------------------------------------------------------------------
         IF (Y90WAF(PVROW,'P')) THEN
            DO 100 I = M, 1, -1
               CALL F06EGF(N,A(I,1),IA,A(IPVROW(I),1),IA)
  100       CONTINUE
         END IF
C
         IF (Y90WAF(PVCOL,'P')) THEN
            DO 120 I = N, 1, -1
               CALL F06EGF(N,A(1,I),1,A(1,IPVCOL(I)),1)
  120       CONTINUE
         END IF
C
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90RCF
C
C-----------------------------------------------------------------------
      RETURN
      END
