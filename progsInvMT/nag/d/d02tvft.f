      SUBROUTINE D02TVF(NEQ,M,NLBC,NRBC,NCOL,TOLS,MXMESH,NMESH,MESH,
     *                  IPMESH,RWORK,LRWORK,IWORK,LIWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02TVF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LIWORK, LRWORK, MXMESH, NCOL, NEQ, NLBC,
     *                  NMESH, NRBC
C     .. Array Arguments ..
      DOUBLE PRECISION  MESH(MXMESH), RWORK(LRWORK), TOLS(NEQ)
      INTEGER           IPMESH(MXMESH), IWORK(LIWORK), M(NEQ)
C     .. Local Scalars ..
      INTEGER           I, IER, IIP, IRP, JP, KD, LIREQ, LRREQ, MATCH,
     *                  MAXORD, MSTAR, MXFIXP, NFIXP, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IWORK(1) = -2
      IER = 0
      NREC = 0
      IF (NEQ.LT.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,I6,A)') ' ** You have set NEQ = ', NEQ,
     *     ' which is less than 1.'
      ELSE IF (NMESH.LT.6) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,I6,A)') ' ** You have set NMESH = ', NMESH,
     *     ' which is less than 6.'
      ELSE IF (MXMESH.LT.2*NMESH-1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,I6,A,I6,A)') ' ** You have set MXMESH = ',
     *     MXMESH, ' which is less than 2*NMESH-1 = ', 2*NMESH - 1, '.'
      END IF
      IF (IER.EQ.1) GO TO 140
C
      MSTAR = 0
      MAXORD = 0
      DO 20 I = 1, NEQ
         IF (M(I).LT.1 .OR. M(I).GT.4) THEN
            IER = 1
            NREC = 1
            WRITE (REC,FMT='(A,I3,A,I6,A)')
     *        ' ** You have set component ', I, ' of M to be ', M(I),
     *        ' which is not 1,2,3 or 4.'
         ELSE IF (TOLS(I).LT.100.0D0*X02AJF()) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(a,i3,a,e12.4/a,e12.4,a)')
     *        ' ** You have set component ', I, ' of TOLS to be ',
     *        TOLS(I), ' ** which is less than the permitted minimum',
     *        100.0D0*X02AJF(), '.'
         ELSE
            MSTAR = MSTAR + M(I)
            MAXORD = MAX(MAXORD,M(I))
         END IF
         IF (IER.EQ.1) GO TO 140
   20 CONTINUE
C
      IF (NCOL.LT.MAXORD .OR. NCOL.GT.7) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,i3,a/a,i1,a,a)')
     *     ' ** You have set NCOL to be ', NCOL, ' which is invalid.',
     *     ' ** It must be .ge. max(M(i),i=1,NEQ), = ', MAXORD, ',',
     *     ' and .le. 7.'
      ELSE IF (NLBC.LT.1 .OR. NRBC.LT.1) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,i6,a,i6,a/a)')
     *     ' ** You have set NLBC and NRBC to be', NLBC, ' and', NRBC,
     *     '.', ' ** They must both be .ge. 1.'
      ELSE IF (NLBC+NRBC.NE.MSTAR) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,i6,a/a,i6,a)')
     *     ' ** You have set NLBC and NRBC such that their sum is',
     *     NLBC + NRBC, '.',
     *     ' ** It must be equal to sum(M(i),i=1,NEQ), which is', MSTAR,
     *     '.'
      END IF
      IF (IER.EQ.1) GO TO 140
C
      DO 40 I = 1, NMESH - 1
         IF (MESH(I).GE.MESH(I+1)) THEN
            IER = 1
            NREC = 1
            WRITE (REC,FMT='(a,a)')
     *        ' ** You have not set the entries of MESH to be',
     *        ' strictly increasing.'
            GO TO 140
         END IF
   40 CONTINUE
C
      IF (IPMESH(1).NE.1 .OR. IPMESH(NMESH).NE.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(a)')
     *     ' ** You have not set IPMESH(1) and IPMESH(NMESH) to be 1.'
         GO TO 140
      END IF
C
      NFIXP = 0
      DO 60 I = 2, NMESH - 1
         IF (IPMESH(I).EQ.1) THEN
            NFIXP = NFIXP + 1
         ELSE IF (IPMESH(I).NE.2) THEN
            IER = 1
            NREC = 1
            WRITE (REC,FMT='(a,a)')
     *        ' ** You have not set the entries of IPMESH to be',
     *        ' 1 or 2.'
            GO TO 140
         END IF
   60 CONTINUE
C
      KD = NEQ*NCOL
      LRREQ = 50 + NEQ*(6+MAXORD*(1+NEQ+MAX(NLBC,NRBC))) -
     *        KD*(KD+6+MSTAR) - MSTAR*(MSTAR-2) +
     *        MXMESH*(2*MSTAR**2+9*MSTAR+6+KD*(KD+MSTAR+6)) + MXMESH/2
      LIREQ = 3*NEQ + 23 - KD + MXMESH*(MSTAR+KD+4)
C      LIREQ = 6 + MSTAR + MXMESH*(MSTAR+KD+4) + 3*NEQ + 21
      IF (LRREQ.GT.LRWORK) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,a/a,i8,a,i8,a)')
     *     ' ** You have not supplied sufficient real workspace given',
     *     ' the input parameters.', ' ** LRWORK was given as ', LRWORK,
     *     ' but should be ', LRREQ, '.'
      ELSE IF (LIREQ.GT.LIWORK) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(a,a/a,i6,a,i6,a)')
     *     ' ** You have not supplied sufficient integer workspace',
     *     ' given the input', ' ** parameters. LIWORK was given as ',
     *     LIWORK, ' but should be ', LIREQ, '.'
      END IF
      IF (IER.EQ.1) GO TO 140
C
      IWORK(1) = 2
      IWORK(2) = NEQ
      IWORK(3) = NCOL
      IWORK(4) = NMESH - 1
      IWORK(5) = MXMESH - 1
      IWORK(6) = MAXORD
      IWORK(7) = MSTAR
      IWORK(8) = NFIXP
      IWORK(9) = NLBC
      IWORK(10) = 1
C      iwork(20) = mxfixp
      MXFIXP = MXMESH/2 + 1
      IWORK(20) = MXFIXP
C
      IWORK(1) = IER + 1
      IWORK(12) = -1
      IWORK(13) = -1
      IWORK(14) = -1
      IWORK(15) = -1
C  Ensures no debugging messages printed - NAG use only
      IWORK(17) = 1
C
      IIP = 20
      JP = 1
      DO 80 I = 1, NEQ
         IWORK(IIP+I) = M(I)
         IWORK(IIP+NEQ+I) = I
         IWORK(IIP+2*NEQ+I) = JP
         JP = JP + M(I)
   80 CONTINUE
      IIP = IIP + 3*NEQ
C
C      DO 140 I = 1, NFIXP
C         RWORK(I) = FIXP(I)
C  140 CONTINUE

      IRP = MXFIXP + 1
      DO 100 I = 1, NEQ
         RWORK(IRP+I) = TOLS(I)
  100 CONTINUE
      IRP = IRP + NEQ + MXMESH*(2+MSTAR) + 49 + KD*(MXMESH-1)
C
      MATCH = 1
      RWORK(IRP+1) = MESH(1)
      DO 120 I = 2, NMESH - 1
         RWORK(IRP+I) = MESH(I)
         IF (IPMESH(I).EQ.1) THEN
            RWORK(MATCH+1) = MESH(I)
            MATCH = MATCH + 1
         END IF
  120 CONTINUE
      RWORK(IRP+NMESH) = MESH(NMESH)
      IRP = IRP + MXMESH + 1
C
  140 CONTINUE
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
      END
