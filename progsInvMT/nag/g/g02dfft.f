      SUBROUTINE G02DFF(IP,Q,LDQ,INDX,RSS,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     REMOVES THE INDX TH X VARIABLE FROM THE REGRESSION.
C     ONLY R AND Q'Y ARE UPDATED.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02DFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS
      INTEGER           IFAIL, INDX, IP, LDQ
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,IP+1), WK(2*IP)
C     .. Local Scalars ..
      INTEGER           I, IERROR, K, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FBF, F06QRF, F06QXF, DCOPY
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (IP.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) IP
      ELSE IF (LDQ.LT.IP) THEN
         WRITE (P01REC(1),FMT=99998) LDQ, IP
      ELSE IF (INDX.LT.1 .OR. INDX.GT.IP) THEN
         WRITE (P01REC(1),FMT=99997) INDX
      ELSE IF (RSS.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99995) RSS
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         K = IP - INDX
         IF (K.NE.0) THEN
C
C     MOVE INDX COLUMN OF R TO END AND SHUFFLE COLUMNS
C     ALONG
C
            CALL DCOPY(INDX,Q(1,INDX+1),1,WK(IP+1),1)
            DO 20 I = INDX + 1, IP
               IF (Q(I,I+1).EQ.0.0D0) THEN
                  GO TO 40
C
               ELSE
                  WK(I-1) = Q(I,I+1)
                  CALL DCOPY(I-1,Q(1,I+1),1,Q(1,I),1)
                  Q(I,I+1) = 0.0D0
               END IF
   20       CONTINUE
C
C     PLACE INDX COLUMN AT THE END OF R
C
            CALL DCOPY(INDX,WK(IP+1),1,Q(1,IP+1),1)
            CALL F06FBF(IP-INDX,0.0D0,Q(INDX+1,IP+1),1)
C
C     RESTORE R TO UPPER TRIANGULAR FORM
C
            CALL F06QRF('L',IP,INDX,IP,WK(IP+1),WK,Q(1,2),LDQ)
            CALL F06QXF('L','V','F',IP,1,INDX,IP,WK(IP+1),WK,Q,IP)
            GO TO 60
C
   40       IERROR = 2
            WRITE (P01REC(1),FMT=99996) I, I + 1
            GO TO 80
C
         END IF
   60    RSS = Q(IP,1)*Q(IP,1) + RSS
      END IF
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, IP.lt.1 : IP = ',I16)
99998 FORMAT (' ** On entry, LDQ.lt.IP : LDQ = ',I16,' IP = ',I16)
99997 FORMAT (' ** On entry, INDX.lt.1 or INDX.gt.IP : INDX = ',I16)
99996 FORMAT (' ** On entry, Q(',I16,',',I16,').eq.0.0')
99995 FORMAT (' ** On entry, RSS.lt.0.0 : RSS = ',D13.5)
      END
