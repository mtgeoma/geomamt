      SUBROUTINE D02TXF(MXMESH,NMESH,MESH,IPMESH,RWORK,IWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02TXF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, MXMESH, NMESH
C     .. Array Arguments ..
      DOUBLE PRECISION  MESH(MXMESH), RWORK(*)
      INTEGER           IPMESH(MXMESH), IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, XR
      INTEGER           COUNT, I, IC, IER, IMSH, KD, MATCH, MSTAR,
     *                  MXFIXP, NCOL, NEQ, NFIXP, NMAX, NOLD, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
      NEQ = IWORK(2)
      NCOL = IWORK(3)
      KD = NCOL*NEQ
      NOLD = IWORK(4)
      MSTAR = IWORK(7)
      MXFIXP = IWORK(20)
      NMAX = IWORK(5)
      IMSH = 1 + NEQ + MXFIXP
      IF (IWORK(16).LT.1 .OR. IWORK(16).GT.6) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(a)')
     *     ' ** The solver routine does not appear to have been called.'
      ELSE IF (IWORK(16).GT.1 .AND. IWORK(16).LE.6) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a/a)')
     * ' ** The solver routine did not produce any results suitable for'
     *     , ' ** remeshing.'
      ELSE IF (MXMESH.NE.NMAX+1) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,i6,a,a/a,i6,a)')
     *     ' ** The value of MXMESH = ', MXMESH, ', which is not that',
     *     ' supplied to', ' ** the setup routine = ', NMAX + 1, '.'
         GO TO 100
      ELSE IF (NMESH.LT.2) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,I6,A,I6,a)') ' ** You have set NMESH = ',
     *     NMESH, ' which is less than 2.'
      ELSE IF (NMESH.GT.(MXMESH-1)/2) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,I6,A,I6,a)') ' ** You have set NMESH = ',
     *     NMESH, ' which is greater than the maximum,', (MXMESH-1)/2,
     *     '.'
      ELSE IF (MESH(1).NE.RWORK(IMSH+1)) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,e12.4,a,a/a,e12.4,a)') ' ** MESH(1), ',
     *     MESH(1), ', does not coincide with the ', 'left hand end of',
     *     ' ** the range previosuly specified, ', RWORK(IMSH+1), '.'
      ELSE IF (IPMESH(1).NE.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(a,i6,a)') ' ** IPMESH(1) = ', IPMESH(1),
     *     ', which is not 1.'
      ELSE
         COUNT = 0
         NFIXP = 0
         DO 20 I = 1, MXMESH
            IC = IPMESH(I)
            IF (IC.EQ.1) THEN
               COUNT = COUNT + 1
               NFIXP = NFIXP + 1
               X = MESH(I)
            ELSE IF (IC.EQ.2) THEN
               COUNT = COUNT + 1
               X = MESH(I)
            ELSE IF (IC.NE.3 .AND. IC.NE.-1) THEN
               IER = 1
               NREC = 1
               WRITE (REC,FMT='(a)')
     * ' ** You have set an element of IPMESH not equal to -1,1,2 or 3.'
               GO TO 100
            END IF
            IF (COUNT.EQ.NMESH .OR. IC.EQ.-1) GO TO 40
   20    CONTINUE
C
   40    CONTINUE
         IF (IC.EQ.-1) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(a/a)')
     *    ' ** An element of IPMESH was set to -1 before NMESH elements'
     *        , ' ** containing 1 or 2 were detected.'
         ELSE IF (COUNT.NE.NMESH) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(a,i6,a/a,i6,a)') ' ** Expected ', NMESH,
     *        ' elements of IPMESH to be 1 or 2, but',
     *        ' ** only found ', COUNT, '.'
         ELSE IF (IC.NE.1) THEN
            IER = 1
            NREC = 3
            WRITE (REC,FMT='(a/a/a,i6,a)')
     *        ' ** You have set the element of IPMESH correpsonding to',
     *        ' ** the last element of MESH to be included in the new',
     *        ' ** new mesh as ', IC, ', which is not 1.'
         ELSE IF (X.NE.RWORK(IMSH+NOLD+1)) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(a,e12.4,a/a,e12.4,a)')
     *        ' ** The last point of the new mesh,', X,
     *        ' does not coinicide',
     *  ' ** with the right hand end of the range previously specified,'
     *        , RWORK(IMSH+NOLD+1), '.'
         ELSE
C
C next line signifies mesh change
            IWORK(12) = 1
            IWORK(10) = 4
            NFIXP = NFIXP - 2
            IWORK(8) = NFIXP
            IWORK(4) = NMESH - 1
C            LENCP = (NOLD+1) + (NOLD+1)*MSTAR + NEQ*NCOL*NOLD +
C     *              NCOL*NCOL
C            IMSH = 1 + MXFIXP + NEQ
C            IOLD = IMSH + LENCP + 1
C            INEW = IOLD + NMESH
C            DO 60 I = 1, LENCP
C               RWORK(INEW-I) = RWORK(IOLD-I)
C   60       CONTINUE
C
            IMSH = 1 + MXFIXP + NEQ + NMAX + 1 + MSTAR*(NMAX+1) +
     *             KD*NMAX + 49
            X = MESH(1)
            RWORK(IMSH+1) = X
            MATCH = 0
            COUNT = 1

            DO 60 I = 2, MXMESH
               IC = IPMESH(I)
               IF (IC.EQ.1 .OR. IC.EQ.2) THEN
                  COUNT = COUNT + 1
                  XR = MESH(I)
                  IF (XR.LE.X) THEN
                     IER = 1
                     NREC = 1
                     WRITE (REC,FMT='(a,a)')
     *         ' ** You have not set the entries of MESH to be strictly'
     *                 , ' increasing.'
                     GO TO 100
                  ELSE
                     IF (MATCH.LT.NFIXP) THEN
                        IF (IC.EQ.1) THEN
                           MATCH = MATCH + 1
                           RWORK(MATCH+1) = XR
                        END IF
                     END IF
                     RWORK(IMSH+COUNT) = XR
                     X = XR
                  END IF
               END IF
               IF (COUNT.EQ.NMESH) GO TO 80
   60       CONTINUE
   80       CONTINUE
         END IF
      END IF
C
  100 CONTINUE
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
      END
