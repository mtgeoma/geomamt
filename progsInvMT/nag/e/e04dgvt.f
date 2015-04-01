      LOGICAL FUNCTION E04DGV(OBJF,OLDF,ALFA,PNORM,RTFTOL,CUBERT,FTOL,
     *                        XNORM,GNORM,EPSAF)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     E04DGV returns the value true if the CG method has converged.
C
C     Refer to Gill et al, practical optimization pages 306-307
C     for termination criteria.
C
C     FTOL here corresponds to the tau of formulas 1 - 4 there.
C
C     -- Written on 4-June-1986.
C     Sven Hammarling and Janet Welding, NAG Central Office.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        ALFA, CUBERT, EPSAF, FTOL, GNORM, OBJF,
     *                        OLDF, PNORM, RTFTOL, XNORM
C     .. Arrays in Common ..
      DOUBLE PRECISION        WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION        FTEST, THETA
C     .. Intrinsic Functions ..
      INTRINSIC               ABS
C     .. Common blocks ..
      COMMON                  /AX02ZA/WMACH
C     .. Save statement ..
      SAVE                    /AX02ZA/
C     .. Executable Statements ..
      IF (GNORM.LT.EPSAF) THEN
         E04DGV = .TRUE.
      ELSE
         FTEST = 1 + ABS(OBJF)
         THETA = FTOL*FTEST
         IF ((OLDF-OBJF).GE.THETA) THEN
            E04DGV = .FALSE.
         ELSE IF ((ALFA*PNORM).GE.(RTFTOL*(1+XNORM))) THEN
            E04DGV = .FALSE.
         ELSE IF (GNORM.GT.(CUBERT*FTEST)) THEN
            E04DGV = .FALSE.
         ELSE
            E04DGV = .TRUE.
         END IF
      END IF
      RETURN
C
C     End of E04DGV (CNVERG).
C
      END
