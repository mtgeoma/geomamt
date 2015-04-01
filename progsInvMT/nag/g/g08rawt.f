      SUBROUTINE G08RAW(PAREST,NPEST,PARVAR,NPVAR,IP,WA,NWA,WB,IERROR)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATES PARAMETER ESTIMATES AND CHI-SQUARED STATISTICS
C     FOR THE LINEAR MODEL USING A LIKELIHOOD BASED ON RANKS.
C
C     PETTITT, A.N.P. INFERENCE FOR THE LINEAR MODEL USING A
C                     LIKELIHOOD BASED ON RANKS.
C
C     ARGUMENTS :
C              PAREST - REAL ARRAY OF DIMENSION AT LEAST NPEST.
C                       ON ENTRY FIRST IP ELEMENTS CONTAIN
C                       SCORE STATISTIC. ON EXIT ALSO CONTAINS
C                       PARAMETER ESTIMATES, CHI-SQUARED STATISTIC,
C                       S.E.'S OF PARAMETER ESTIMATES AND Z-STATISTICS.
C               NPEST - INTEGER SPECIFYING LENGTH OF PAREST.
C                           (NPEST .GE. 4*IP+1)
C              PARVAR - REAL ARRAY OF DIMENSION AT LEAST (IP+1,IP).
C                       ON ENTRY CONTAINS VAR-COVARIANCE MATRIX OF
C                       SCORE STATISTIC IN UPPER TRIANGLE. ON EXIT
C                       CONTAINS VAR-COVARIANCE MATRIX OF PARAMETER
C                       ESTIMATES IN LOWER TRIANGLE.
C               NPVAR - INTEGER SPECIFYING FIRST DIMENSION OF PARVAR
C                       AS DEFINED IN CALLING (SUB)PROGRAM.
C                  IP - INTEGER SPECIFYING NUMBER OF PARAMETERS FITTED.
C                  WA - REAL ARRAY OF DIMENSION AT LEAST (IP,IP).
C                 NWA - FIRST DIMENSION OF WA AS DEFINED IN CALLING
C                       CALLING (SUB)PROGRAM. NWA .GE. IP.
C                       USED AS WORKSPACE.
C                  WB - REAL ARRAY OF DIMENSION AT LEAST (IP). USED
C                       AS WORKSPACE.
C              IERROR - IF IERROR IS NOT EQUAL TO 0 ON EXIT, THE
C                        MATRIX PARVAR CANNOT BE INVERTED.
C
C
C     FIRST INVERT MATRIX
C
C     .. Scalar Arguments ..
      INTEGER           IERROR, IP, NPEST, NPVAR, NWA
C     .. Array Arguments ..
      DOUBLE PRECISION  PAREST(NPEST), PARVAR(NPVAR,IP), WA(NWA,IP),
     *                  WB(IP)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, TERM
      INTEGER           I, IFAIL, J, K, L, MPOS
C     .. External Subroutines ..
      EXTERNAL          F01ABF
C     .. Executable Statements ..
      IFAIL = 1
      CALL F01ABF(PARVAR,NPVAR,IP,WA,NWA,WB,IFAIL)
      IF (IFAIL.NE.0) THEN
         IERROR = 4
         RETURN
      END IF
C
C     CALCULATE PARAMETER ESTIMATES
C
      DO 40 I = 1, IP
         DO 20 J = 1, IP
            IF (I.GE.J) THEN
               K = I + 1
               L = J
            ELSE
               K = J + 1
               L = I
            END IF
            PAREST(IP+I) = PAREST(IP+I) + PARVAR(K,L)*PAREST(J)
   20    CONTINUE
   40 CONTINUE
C
C     CALCULATE CHI-SQUARED STATISTIC
C
      SUM = 0.0D0
      DO 80 I = 1, IP
         TERM = 0.0D0
         DO 60 J = 1, IP
            IF (I.LE.J) THEN
               K = I
               L = J
            ELSE
               K = J
               L = I
            END IF
            TERM = TERM + PARVAR(K,L)*PAREST(IP+J)
   60    CONTINUE
         SUM = SUM + TERM*PAREST(IP+I)
   80 CONTINUE
      MPOS = 2*IP + 1
      PAREST(MPOS) = SUM
C
C     CALCULATE STANDARD ERRORS AND Z-STATISTICS
C
      DO 100 I = 1, IP
         PAREST(MPOS+I) = PARVAR(I+1,I)**0.5D0
         PAREST(MPOS+IP+I) = PAREST(MPOS-IP-1+I)/PAREST(MPOS+I)
  100 CONTINUE
      RETURN
      END
