      SUBROUTINE F01BSZ(N,A,LICN,IVECT,JVECT,NZ,ICN,LENR,LENRL,LENOFF,
     *                  IP,IQ,IW1,IW,W1,ABORT,IDISP,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7E REVISED IER-202 (JUL 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DERIVED FROM HARWELL LIBRARY ROUTINE MA28D
C
C     SORTS THE USER'S MATRIX INTO THE STRUCTURE OF THE DECOMPOSED
C     FORM AND CHECKS FOR THE PRESENCE OF DUPLICATE ENTRIES OR
C     NON-ZEROS LYING OUTSIDE THE SPARSITY PATTERN OF THE
C     DECOMPOSITION.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  W1
      INTEGER           IFAIL, LICN, N, NZ
      LOGICAL           ABORT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LICN)
      INTEGER           ICN(LICN), IDISP(2), IP(N), IQ(N), IVECT(NZ),
     *                  IW(N,2), IW1(N,3), JVECT(NZ), LENOFF(N),
     *                  LENR(N), LENRL(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, ZERO
      INTEGER           I, IBLOCK, IDISP2, IDUMMY, II, INEW, IOLD,
     *                  ISAVE, J1, J2, JCOMP, JDUMMY, JJ, JNEW, JOLD,
     *                  MIDPT, NADV, NERR
      LOGICAL           BLOCKL
C     .. Local Arrays ..
      CHARACTER*65      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MOD
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 0
C     NERR IS THE UNIT NUMBER FOR ERROR MESSAGES
C     NADV IS THE UNIT NUMBER FOR DUPLICATE ELEMENT WARNING MESSAGES
      CALL X04AAF(0,NERR)
      CALL X04ABF(0,NADV)
      BLOCKL = LENOFF(1) .GE. 0
      IBLOCK = 1
C     IW1(I,3)  IS SET TO THE BLOCK IN WHICH ROW I LIES AND THE
C     INVERSE PERMUTATIONS TO IP AND IQ ARE SET IN IW1(.,1) AND
C     IW1(.,2) RESP.
C     POINTERS TO BEGINNING OF THE PART OF ROW I IN DIAGONAL AND
C     OFF-DIAGONAL BLOCKS ARE SET IN IW(I,2) AND IW(I,1) RESP.
      IW(1,1) = 1
      IW(1,2) = IDISP(1)
      DO 20 I = 1, N
         IW1(I,3) = IBLOCK
         IF (IP(I).LT.0) IBLOCK = IBLOCK + 1
         II = ABS(IP(I))
         IW1(II,1) = I
         JJ = IQ(I)
         JJ = ABS(JJ)
         IW1(JJ,2) = I
         IF (I.EQ.1) GO TO 20
         IF (BLOCKL) IW(I,1) = IW(I-1,1) + LENOFF(I-1)
         IW(I,2) = IW(I-1,2) + LENR(I-1)
   20 CONTINUE
C     PLACE EACH NON-ZERO IN TURN INTO ITS CORRECT LOCATION
C     IN THE A/ICN ARRAY.
      IDISP2 = IDISP(2)
      DO 340 I = 1, NZ
         IF (I.GT.IDISP2) GO TO 40
         IF (ICN(I).LT.0) GO TO 340
   40    IOLD = IVECT(I)
         JOLD = JVECT(I)
         AA = A(I)
C        THIS IS A DUMMY LOOP FOR FOLLOWING A CHAIN OF INTERCHANGES.
C        IT WILL BE EXECUTED NZ TIMES IN TOTAL.
         DO 280 IDUMMY = 1, NZ
C           PERFORM SOME VALIDITY CHECKS ON IOLD AND JOLD.
            IF (IOLD.LE.N .AND. IOLD.GT.0 .AND. JOLD.LE.N .AND. JOLD.GT.
     *          0) GO TO 60
            IFAIL = 4
C           ** CODE FOR OUTPUT OF ERROR MESSAGE **
            IF (MOD(ISAVE/10,10).NE.0) THEN
               WRITE (REC,FMT=99999) I, A(I), IOLD, JOLD
               CALL X04BAF(NERR,REC(1))
               CALL X04BAF(NERR,REC(2))
               CALL X04BAF(NERR,REC(3))
            END IF
C           ** END OF CODE FOR OUTPUT OF ERROR MESSAGE **
            GO TO 360
   60       INEW = IW1(IOLD,1)
            JNEW = IW1(JOLD,2)
C           ARE WE IN A VALID BLOCK AND IS IT DIAGONAL OR
C           OFF-DIAGONAL...
            IF (IW1(INEW,3)-IW1(JNEW,3)) 80, 120, 100
   80       IFAIL = 5
C           ** CODE FOR OUTPUT OF ERROR MESSAGE **
            IF (MOD(ISAVE/10,10).NE.0) THEN
               WRITE (REC,FMT=99998) IOLD, JOLD
               CALL X04BAF(NERR,REC(1))
               CALL X04BAF(NERR,REC(2))
            END IF
C           ** END OF CODE FOR OUTPUT OF ERROR MESSAGE **
            GO TO 360
  100       J1 = IW(INEW,1)
            J2 = J1 + LENOFF(INEW) - 1
            GO TO 220
C           ELEMENT IS IN DIAGONAL BLOCK.
  120       J1 = IW(INEW,2)
            IF (INEW.GT.JNEW) GO TO 140
            J2 = J1 + LENR(INEW) - 1
            J1 = J1 + LENRL(INEW)
            GO TO 220
  140       J2 = J1 + LENRL(INEW)
C           BINARY SEARCH OF ORDERED LIST  .. ELEMENT IN L PART OF ROW.
            DO 200 JDUMMY = 1, N
               MIDPT = (J1+J2)/2
               JCOMP = ABS(ICN(MIDPT))
               IF (JNEW-JCOMP) 160, 260, 180
  160          J2 = MIDPT
               GO TO 200
  180          J1 = MIDPT
  200       CONTINUE
            IFAIL = 5
C           ** CODE FOR OUTPUT OF ERROR MESSAGE **
            IF (MOD(ISAVE/10,10).NE.0) THEN
               WRITE (REC,FMT=99997) IOLD, JOLD
               CALL X04BAF(NERR,REC(1))
               CALL X04BAF(NERR,REC(2))
            END IF
C           ** END OF CODE FOR OUTPUT OF ERROR MESSAGE **
            GO TO 360
C           LINEAR SEARCH ... ELEMENT IN U PART OF ROW OR OFF-DIAGONAL
C           BLOCKS.
  220       DO 240 MIDPT = J1, J2
               IF (ABS(ICN(MIDPT)).EQ.JNEW) GO TO 260
  240       CONTINUE
            IFAIL = 5
C           ** CODE FOR OUTPUT OF ERROR MESSAGE **
            IF (MOD(ISAVE/10,10).NE.0) THEN
               WRITE (REC,FMT=99997) IOLD, JOLD
               CALL X04BAF(NERR,REC(1))
               CALL X04BAF(NERR,REC(2))
            END IF
C           ** END OF CODE FOR OUTPUT OF ERROR MESSAGE **
            GO TO 360
C           EQUIVALENT ELEMENT OF ICN IS IN POSITION MIDPT.
  260       IF (ICN(MIDPT).LT.0) GO TO 320
            IF (MIDPT.GT.NZ .OR. MIDPT.LE.I) GO TO 300
            W1 = A(MIDPT)
            A(MIDPT) = AA
            AA = W1
            IOLD = IVECT(MIDPT)
            JOLD = JVECT(MIDPT)
            ICN(MIDPT) = -ICN(MIDPT)
  280    CONTINUE
  300    A(MIDPT) = AA
         ICN(MIDPT) = -ICN(MIDPT)
         GO TO 340
C        DUPLICATE ELEMENT FOUND
  320    A(MIDPT) = A(MIDPT) + AA
C        ** CODE FOR OUTPUT OF WARNING MESSAGE ************************
         IF (MOD(ISAVE/100,100).NE.0) THEN
            WRITE (REC,FMT=99996) IOLD, JOLD, AA
            CALL X04BAF(NADV,REC(1))
            CALL X04BAF(NADV,REC(2))
            CALL X04BAF(NADV,REC(3))
         END IF
C        ** END OF CODE FOR OUTPUT OF WARNING MESSAGE *****************
         IF (ABORT) IFAIL = 8
  340 CONTINUE
C     RESET ICN ARRAY  AND ZERO ELEMENTS IN L/U BUT NOT IN A.
C     ALSO CALCULATE MAXIMUM ELEMENT OF A.
  360 W1 = ZERO
      DO 400 I = 1, IDISP2
         IF (ICN(I).LT.0) GO TO 380
         A(I) = ZERO
         GO TO 400
  380    ICN(I) = -ICN(I)
         W1 = MAX(W1,ABS(A(I)))
  400 CONTINUE
      RETURN
C
99999 FORMAT (/' ON ENTRY IRN(I) OR ICN(I) IS OUT OF RANGE - I = ',I8,
     *  /'  A(I) = ',1P,D14.6,'  IRN(I) = ',I8,'  ICN(I) = ',I8)
99998 FORMAT (/' NON-ZERO ELEMENT (',I6,',',I6,') IN ZERO OFF-DIAGONAL',
     *  ' BLOCK')
99997 FORMAT (/' NON-ZERO ELEMENT (',I6,',',I6,') WAS NOT IN L/U PATTE',
     *  'RN')
99996 FORMAT (/' F01BSF FOUND DUPLICATE ELEMENT WITH INDICES ',I6,', ',
     *  I6,/'  VALUE = ',1P,D14.6)
      END
