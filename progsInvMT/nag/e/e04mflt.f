      SUBROUTINE E04MFL(MSGLVL,N,NRZ,NZ,ZEROLM,NOTOPT,NUMINF,TRUSML,
     *                  SMLLST,JSMLST,TINYST,JTINY,GQ)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1558 (JUN 1995).
C
C     ******************************************************************
C     E04MFL  updates JSMLST and SMLLST when there are artificial
C     constraints.
C
C     On input,  JSMLST  is the index of the minimum of the set of
C     adjusted multipliers.
C     On output, a negative JSMLST defines the index in Q'g of the
C     artificial constraint to be deleted.
C
C     Original version written 17-Jan-1988.
C     This version of E04MFL dated  23-Jul-1991.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SMLLST, TINYST, TRUSML, ZEROLM
      INTEGER           JSMLST, JTINY, MSGLVL, N, NOTOPT, NRZ, NUMINF,
     *                  NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  GQ(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  RLAM
      INTEGER           J, K, KK, LENGTH
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      DO 20 J = NRZ + 1, NZ
         RLAM = -ABS(GQ(J))
C
         IF (RLAM.LT.ZEROLM) THEN
            IF (NUMINF.EQ.0) NOTOPT = NOTOPT + 1
C
            IF (RLAM.LT.SMLLST) THEN
               TRUSML = GQ(J)
               SMLLST = RLAM
               JSMLST = -J
            END IF
C
         ELSE IF (RLAM.LT.TINYST) THEN
            TINYST = RLAM
            JTINY = -J
         END IF
   20 CONTINUE
C
      IF (MSGLVL.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999)
            CALL X04BAY(ISUMM,2,REC)
            LENGTH = NZ - NRZ
            DO 40 K = 1, LENGTH, 4
               WRITE (REC,FMT=99998) (GQ(KK),KK=K,MIN(K+3,LENGTH))
               CALL X04BAF(ISUMM,REC(1))
   40       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of  E04MFL.  (CMMUL2)
C
99999 FORMAT (/' Multipliers for the artificial constraints        ')
99998 FORMAT (4(5X,1P,D11.2))
      END
