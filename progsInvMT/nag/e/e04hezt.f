      SUBROUTINE E04HEZ(M,N,FVEC,FJAC,LJ,G)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     THIS SUBROUTINE FORMS THE GRADIENT VECTOR OF THE SUM OF
C     SQUARES FOR LEAST-SQUARES PROBLEMS.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN
C     AND NICHOLAS I. M. GOULD.
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           LJ, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), G(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, J
C     .. Executable Statements ..
      DO 40 I = 1, N
         SUM = 0.0D+0
         DO 20 J = 1, M
            SUM = SUM + FJAC(J,I)*FVEC(J)
   20    CONTINUE
         G(I) = SUM + SUM
   40 CONTINUE
      RETURN
C
C     END OF E04HEZ   (LSQGRD)
C
      END
