      SUBROUTINE E01SEF(M,X,Y,F,RNW,RNQ,NW,NQ,FNODES,MINNQ,WRK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14C REVISED. IER-876 (NOV 1990).
C
C     This subroutine serves to construct a Shepard's method type of
C     surface through a set of scattered data points, using a least
C     squares quadratic nodal function through each point.
C
C     Routine created - December 1986
C     Author          - Richard Franke
C                       Naval Postgraduate School Monterey,
C                       California  93940
C                       Adapted for Nag by H.Scullion (Leic Univ.)
C                       and I. Gladwell (Nag Ltd.)
C
C     Input Parameters:
C
C            M  -  The number of data points.
C
C        X,Y,F  -  The data points, (X(I),Y(I),F(I),I=1,M)
C
C          RNW  -  The radius for the weight functions.
C
C          RNQ  -  The radius for the nodal functions. If either of
C                  RNW or RNQ is .le. zero on entry, then their values
C                  will be computed, using the values of NW and NQ,
C                  and will be returned on exit.
C
C           NW  -  The approximate number of neighbouring nodes to
C                  affect each weight function.
C
C           NQ  -  The approximate number of neighbouring nodes to
C                  affect each nodal function. If either of
C                  NW or NQ is .le. zero on entry, then a default
C                  value will be used.
C
C     Output Parameters:
C
C       FNODES  -  Real array of dimension at least (5*M).
C                  This array is used to store the coefficients for the
C                  nodal functions.
C
C          WRK  -  Real work array of dimension at least (6*M),
C                  used in E01SEZ.
C
C     On exit, if IFAIL = 0, normal return.
C                   = 1, input error. M.lt.3
C                   = 2, input error. RNQ.lt.RNW
C                   = 3, input error. NQ.lt.NW
C                   = 4, input error. The (X(I),Y(I)) points are
C                        not unique for I = 1,2,...,M.
C
C     .. Parameters ..
      INTEGER           NOMNNW, NOMNNQ
      PARAMETER         (NOMNNW=9,NOMNNQ=18)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01SEF')
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RNQ, RNW
      INTEGER           IFAIL, M, MINNQ, NQ, NW
C     .. Array Arguments ..
      DOUBLE PRECISION  F(M), FNODES(5*M), WRK(6*M), X(M), Y(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DIST
      INTEGER           I, IER, J, NNQ, NNW, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01SEZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT
C     .. Executable Statements ..
      IER = 0
      NREC = 0
      IF (M.LT.3) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT=99999) M
      ELSE
         D = ZERO
         DO 40 I = 2, M
            DO 20 J = 1, I - 1
               DIST = (X(I)-X(J))**2 + (Y(I)-Y(J))**2
               IF (DIST.EQ.ZERO) THEN
                  NREC = 2
                  IER = 4
                  WRITE (REC,FMT=99996) I, J, X(I), Y(I)
                  GO TO 60
               END IF
               D = MAX(DIST,D)
   20       CONTINUE
   40    CONTINUE
         IF (RNQ.LE.ZERO) THEN
            IF (NQ.LE.0) THEN
               NNW = NOMNNW
               NNQ = NOMNNQ
            ELSE
               IF (NQ.LT.NW .OR. NW.LE.0) THEN
                  IER = 3
                  NREC = 1
                  WRITE (REC,FMT=99998) NQ, NW
                  GO TO 60
               ELSE
                  NNW = NW
                  NNQ = NQ
               END IF
            END IF
            RNW = SQRT(NNW*D/(4*M))
            RNQ = RNW*SQRT(DBLE(NNQ)/DBLE(NNW))
         ELSE IF (RNW.LE.ZERO .OR. RNQ.LT.RNW) THEN
            IER = 2
            NREC = 1
            WRITE (REC,FMT=99997) RNQ, RNW
            GO TO 60
         END IF
C
C        Compute the nodal functions.
         MINNQ = M
         CALL E01SEZ(M,X,Y,F,RNQ,FNODES,MINNQ,WRK(1),WRK(5*M+1))
      END IF
   60 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** On entry, M .lt. 3: M =',I16,'.')
99998 FORMAT (' ** On entry, NQ .lt. NW: NQ =',I16,', NW =',I16,'.')
99997 FORMAT (' ** On entry, RNQ .lt. RNW: RNQ =',1P,D13.5,'  RNW =',
     *       D13.5)
99996 FORMAT (1X,'** On entry, data points (X(I),Y(I)) and (X(J),Y(J))',
     *       ' are identical,',/4X,'for I =',I8,', J =',I8,', (X(I),Y(',
     *       'I)) = (',1P,D13.5,',',D13.5,')')
      END
