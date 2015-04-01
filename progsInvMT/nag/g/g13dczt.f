      SUBROUTINE G13DCZ(HESS,N4,PARHLD,DISP,IDISP)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER           IDISP, N4
C     .. Array Arguments ..
      DOUBLE PRECISION  DISP(IDISP,N4), HESS(N4,N4)
      LOGICAL           PARHLD(N4)
C     .. Local Scalars ..
      INTEGER           I, I2, L, L2
C     .. Executable Statements ..
C
C     THIS ROUTINE INSERTS ROWS (AND COLUMNS) OF ZEROS IN THE
C     VARIANCE-COVARIANCE MATRIX TO TAKE ACCOUNT OF PARAMETERS
C     HELD AT THEIR INITIAL VALUES
C
      I2 = 1
      DO 60 I = 1, N4
         IF (PARHLD(I)) THEN
            DO 20 L = 1, N4
               HESS(I,L) = 0.0D0
               HESS(L,I) = 0.0D0
   20       CONTINUE
         ELSE
C
C           SET UP I TH ROW AND COLUMN OF HESS
C
            L2 = 1
            DO 40 L = 1, N4
               IF ( .NOT. PARHLD(L)) THEN
                  HESS(I,L) = DISP(I2,L2)
                  HESS(L,I) = DISP(I2,L2)
                  L2 = L2 + 1
               ELSE
                  HESS(I,L) = 0.0D0
                  HESS(L,I) = 0.0D0
               END IF
C
   40       CONTINUE
C
            I2 = I2 + 1
C
         END IF
C
   60 CONTINUE
C
      RETURN
      END
