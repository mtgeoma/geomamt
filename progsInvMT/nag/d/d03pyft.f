      SUBROUTINE D03PYF(NPDE,U,NBKPTS,XBKPTS,NPOLY,NPTS,XP,INTPTS,ITYPE,
     *                  UOUT,W,NW,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C    -------------------------------------------------------------------
C    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     Routine to provide values of the solution and possibly the
C     first derivative in space and the flux on the mesh
C     XP(INTPTS).
C
C   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   --------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PYF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, INTPTS, ITYPE, NBKPTS, NPDE, NPOLY, NPTS,
     *                  NW
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NPTS*NPDE), UOUT(NPDE,INTPTS,ITYPE), W(NW),
     *                  XBKPTS(NBKPTS), XP(INTPTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, TWOU
      INTEGER           I, I1, I10, I2, IFAIL1, IRESWK, J, K1, NEL,
     *                  NPTL, NPTP, NREC, NVST
C     .. Local Arrays ..
      CHARACTER*80      P01REC(6)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D03PYQ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      IFAIL1 = 0
      TWOU = X02AJF()
C
      IF (NPDE.LE.0) THEN
         NREC = 2
         WRITE (P01REC,FMT=99996) NPDE
         IFAIL1 = 1
         GO TO 100
      END IF
C
      IF (NBKPTS.LE.2) THEN
         NREC = 2
         WRITE (P01REC,FMT=99995) NBKPTS
         IFAIL1 = 1
         GO TO 100
      END IF
C
      IF (NPOLY.LE.0) THEN
         NREC = 2
         WRITE (P01REC,FMT=99994) NPOLY
         IFAIL1 = 1
         GO TO 100
      END IF
C
      IF (INTPTS.LE.0) THEN
         NREC = 2
         WRITE (P01REC,FMT=99993) INTPTS
         IFAIL1 = 1
         GO TO 100
      END IF
C
      IF (NPTS.NE.(NBKPTS-1)*NPOLY+1) THEN
         NREC = 2
         WRITE (P01REC,FMT=99990) NPTS, (NBKPTS-1)*NPOLY + 1
         IFAIL1 = 1
         GO TO 100
      END IF
C
      IF ((XBKPTS(NBKPTS)-XBKPTS(1)).LT.(TWOU*(NBKPTS-1.D0))) THEN
         I1 = 1
         I2 = NBKPTS
         NREC = 5
         IFAIL1 = 1
         WRITE (P01REC,FMT=99992) I1, I2
         GO TO 100
      END IF
C
      DO 20 I = 2, NBKPTS
         IF ((XBKPTS(I)-XBKPTS(I-1)).LT.TWOU) THEN
            I1 = I - 1
            I2 = I
            NREC = 5
            IFAIL1 = 1
            WRITE (P01REC,FMT=99992) I1, XBKPTS(I1), I2, XBKPTS(I2)
            GO TO 100
         END IF
   20 CONTINUE
C
      IRESWK = 1
      NEL = NBKPTS - 1
      NPTL = NPOLY + 1
      NVST = NPDE*NPTS + 1
      NPTP = NPDE*NPTL
      I10 = 3*NPTL*NPTL + NPTL + NPTP + 4*NPDE + 1
C
      IF (ITYPE.LT.1 .OR. ITYPE.GT.2) THEN
         NREC = 2
         WRITE (P01REC,FMT=99999) ITYPE
         IFAIL1 = 1
         GO TO 100
      END IF
C
C     ... Chech the interpolation points XP and if ITYPE = 2 ...
C     ... the points do not conflict with the break-points.  ...
C
      IF ((XP(1).LT.XBKPTS(1)-TWOU) .OR. (XP(INTPTS).GT.XBKPTS(NBKPTS)
     *    +TWOU)) THEN
         NREC = 1
         WRITE (P01REC,FMT=99991)
         IFAIL1 = 3
         GO TO 100
      END IF
C
      DO 40 I = 2, INTPTS
         TEMP = XP(I) - XP(I-1)
         J = I - 1
C
         IF (TEMP.LE.TWOU) THEN
            NREC = 5
            WRITE (P01REC,FMT=99998) J, XP(J), I, XP(I)
            IFAIL1 = 2
            GO TO 100
         END IF
C
   40 CONTINUE
C
      IF (ITYPE.EQ.2 .AND. NBKPTS.GT.2) THEN
         DO 80 I = 1, INTPTS
            DO 60 J = 2, NEL
               TEMP = ABS(XP(I)-XBKPTS(J))
C
               IF (TEMP.LE.TWOU) THEN
                  NREC = 6
                  WRITE (P01REC,FMT=99997) I, XP(I), J, XBKPTS(J)
                  IFAIL1 = 2
                  GO TO 100
               END IF
C
   60       CONTINUE
   80    CONTINUE
      END IF
C
C     ... Call the interpolation routine. ...
C
      K1 = IRESWK + I10 - 1
C
      CALL D03PYQ(INTPTS,XP,UOUT,ITYPE,U,NPTS,NPDE,NEL,NPTL,W(IRESWK),
     *            W(K1),XBKPTS,NBKPTS,IFAIL1)
C
      IF (IFAIL1.EQ.3) THEN
         WRITE (P01REC,FMT=99991)
      END IF
  100 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, ITYPE .ne. 1 or 2, ',/' ** ITYPE =',I16)
99998 FORMAT (' ** On entry, interpolation points are not in ',/' ** s',
     *       'trictly increasing order i.e. point no',I16,/' ** with v',
     *       'alue',D13.5,/' ** is greater than point no',I16,/' ** wi',
     *       'th value',D13.5)
99997 FORMAT (' ** On entry, ITYPE .eq. 2 and at least one interpolati',
     *       'on',/' ** point coincides with a break-point i.e. ',/' *',
     *       '* interpolation point no',I16,/' ** with value ',D13.5,
     *       /' ** is close to break-point',I16,/'** with value',D13.5)
99996 FORMAT (' ** On entry, NPDE .le.  0, ',/' ** NPDE =',I16)
99995 FORMAT (' ** On entry, NBKPTS .le. 2, ',/' ** NBKPTS =',I16)
99994 FORMAT (' ** On entry, NPOLY .le. 0, ',/' ** NPOLY =',I16)
99993 FORMAT (' ** On entry, INTPTS .le.  0, ',/'** INTPTS is : ',I12)
99992 FORMAT (' ** On entry, break-points are not in ',/' ** strictly ',
     *       'increasing order i.e. break-point no',I16,/' ** with val',
     *       'ue',D13.5,/' ** is greater than break-point',I16,/' ** w',
     *       'ith value',D15.5)
99991 FORMAT (' ** Extrapolation is not allowed.')
99990 FORMAT (' ** On entry, NPTS .ne. (NBKPTS-1)*NPOLY + 1',/' ** NPT',
     *       'S =',I16,' it should be',I16)
      END
