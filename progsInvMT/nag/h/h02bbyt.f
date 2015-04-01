      SUBROUTINE H02BBY
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1992.
C     .. Scalars in Common ..
      INTEGER          IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Arrays ..
      CHARACTER*80     REC(2)
C     .. External Subroutines ..
      EXTERNAL         X04BAF, X04BAY
C     .. Common blocks ..
      COMMON           /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      WRITE (REC,FMT=99999)
      CALL X04BAY(NOUT,2,REC)
      WRITE (REC,FMT=99998)
      CALL X04BAF(NOUT,REC(1))
C
      RETURN
C
C
99999 FORMAT (/'  Node  Parent    Obj',7X,'Varbl  Value',6X,'Lower',5X,
     *       'Upper',5X,'Value    Depth')
99998 FORMAT ('   No    Node    Value',6X,'Chosen Before',5X,'Bound',5X,
     *       'Bound',5X,'After')
      END
