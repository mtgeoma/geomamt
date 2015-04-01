      SUBROUTINE C02AFW(ARI,AII,BRI,BII,CR,CI,FAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-868 (NOV 1990).
C
C     Complex Division, (CR,CI) = (AR,AI)/(BR,BI)
C
C     Note: This subroutine is based on the NAG routine F06CLF, but
C           has been modified to avoid complex variables.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AII, ARI, BII, BRI, CI, CR
      LOGICAL           FAIL
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, AR, BI, BR, BIG, DIV, FLMAX, FLMIN, NUMI,
     *                  NUMR, TEMP
      LOGICAL           FIRST
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN
C     .. Save statement ..
      SAVE              BIG, FIRST, FLMAX
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
      AR = ARI
      AI = AII
      BR = BRI
      BI = BII
      IF (AR.EQ.ZERO .AND. AI.EQ.ZERO) THEN
         CR = ZERO
         CI = ZERO
         IF (BR.EQ.ZERO .AND. BI.EQ.ZERO) THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF (FIRST) THEN
            FIRST = .FALSE.
            FLMIN = X02AMF()
            FLMAX = 1/FLMIN
            BIG = FLMAX/2
         END IF
C
         TEMP = MAX(ABS(AR),ABS(AI),ABS(BR),ABS(BI))
         IF (TEMP.GE.BIG) THEN
            AR = AR/2
            AI = AI/2
            BR = BR/2
            BI = BI/2
         END IF
         IF (BR.EQ.ZERO .AND. BI.EQ.ZERO) THEN
            CR = SIGN(FLMAX,AR)
            CI = SIGN(FLMAX,AI)
            FAIL = .TRUE.
         ELSE
            IF (ABS(BR).GE.ABS(BI)) THEN
               TEMP = BI/BR
               DIV = BR + TEMP*BI
               NUMR = AR + TEMP*AI
               NUMI = AI - TEMP*AR
            ELSE
               TEMP = BR/BI
               DIV = BI + TEMP*BR
               NUMR = AI + TEMP*AR
               NUMI = TEMP*AI - AR
            END IF
            IF (ABS(DIV).GE.ONE) THEN
               CR = NUMR/DIV
               CI = NUMI/DIV
               FAIL = .FALSE.
            ELSE
               TEMP = ABS(DIV)*FLMAX
               IF ((ABS(NUMR).LE.TEMP) .AND. (ABS(NUMI).LE.TEMP)) THEN
                  CR = NUMR/DIV
                  CI = NUMI/DIV
                  FAIL = .FALSE.
               ELSE
                  IF (DIV.GE.ZERO) THEN
                     CR = SIGN(FLMAX,NUMR)
                     CI = SIGN(FLMAX,NUMI)
                  ELSE
                     CR = SIGN(FLMAX,-NUMR)
                     CI = SIGN(FLMAX,-NUMI)
                  END IF
                  FAIL = .TRUE.
               END IF
            END IF
         END IF
      END IF
C
      RETURN
      END
