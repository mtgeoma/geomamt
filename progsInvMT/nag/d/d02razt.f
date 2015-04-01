      SUBROUTINE D02RAZ(M,NMAX,N,NUMBEG,NUMMIX,A,B,TOL,X,Y,IY,ABT,FCN,G,
     *                  FCNEP,FCNA,FCNB,JACOBE,JACOBG,JACEPS,JACGEP,A1,
     *                  B1,C1,D1,GAM,A2,B2,WORK,LWORK,IWORK,LIWORK,
     *                  DELEPS,LP,MP,INIT,LIN,IFLAG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-305 (SEP 1981).
C     MARK 9C REVISED. IER-373 (JUN 1982).
C     MARK 10 REVISED. IER-377 (JUN 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-502 (AUG 1986).
C     ***********************************************************
C     AUTHORS
C
C     M. LENTINI  AND  V. PEREYRA  -   OCTOBER 1977           *
C
C     *** REFERENCE ***                               *
C
C     M. LENTINI AND V. PEREYRA ,  AN ADAPTIVE FINITE DIFFE-
C     RENCE SOLVER FOR NONLINEAR TWO-POINT BOUNDARY PROBLEMS
C     WITH MILD BOUNDARY LAYERS. SIAM J. NUMER. ANAL. 14,
C     PP. 91-111 (1977).
C     ***********************************************************
C     FCN, FCNA, FCNB, FCNEP, G, JACEPS, JACGEP, JACOBE, JACOBG
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, DELEPS, TOL
      INTEGER           IFLAG, INIT, IY, LIN, LIWORK, LP, LWORK, M, MP,
     *                  N, NMAX, NUMBEG, NUMMIX
C     .. Array Arguments ..
      DOUBLE PRECISION  A1(M,M), A2(M,2), ABT(M), B1(M), B2(M,2),
     *                  C1(M,M), D1(M,M), GAM(M), WORK(LWORK), X(NMAX),
     *                  Y(IY,NMAX)
      INTEGER           IWORK(LIWORK)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, FCNA, FCNB, FCNEP, G, JACEPS, JACGEP,
     *                  JACOBE, JACOBG
C     .. Local Scalars ..
      INTEGER           I, IM2, IM3, IM4, IM5, IM6, IM7, LWO, M10, M11,
     *                  M12, M14, M15, M16, M22, M23, M4, M5, M6, M7,
     *                  M8, M9, MT, MTNMAX, MTTWO, NADV, NERR, NMAXIW,
     *                  NMAXOK, NMAXW
C     .. Local Arrays ..
      CHARACTER*80      REC(6)
C     .. External Subroutines ..
      EXTERNAL          D02RAY, X04AAF, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IF (LP.NE.0) CALL X04AAF(0,NERR)
      IF (MP.NE.0) CALL X04ABF(0,NADV)
      IFLAG = 0
      IF (M.GE.1 .AND. N.GT.3 .AND. N.LE.NMAX .AND. NUMBEG.GE.0 .AND.
     *    NUMBEG.LT.M .AND. NUMMIX.GE.0 .AND. NUMBEG+NUMMIX.GT.0 .AND.
     *    NUMBEG+NUMMIX.LE.M .AND. NMAX.GE.32) GO TO 40
      IF (LP.NE.0) THEN
         CALL X04BAF(NERR,' ')
         WRITE (REC,FMT=99999) M, N, NMAX, NUMBEG, NUMMIX
         DO 20 I = 1, 6
            CALL X04BAF(NERR,REC(I))
   20    CONTINUE
      END IF
      IFLAG = 1
      GO TO 220
   40 NMAXOK = NMAX
      LWO = LWORK + 2*M*M
      NMAXW = (LWORK-2*M*M-3*M)/(3*M*M+5*M+2)
      IF (LIN.GT.1) NMAXW = (LWORK-2*M*M-3*M)/(3*M*M+6*M+2)
      NMAXIW = (LIWORK-M*M-4*M-2)/(2*M+1)
      IF (LIN.EQ.1 .OR. LIN.EQ.3) NMAXIW = (LIWORK-M)/(2*M+1)
      IF (NMAXW.GE.NMAX .AND. NMAXIW.GE.NMAX) GO TO 160
      NMAXOK = MIN(NMAXW,NMAXIW)
      IF (NMAXOK.GE.N .AND. NMAXOK.GE.32) GO TO 100
      IFLAG = 5
      IF (LP.EQ.0) GO TO 220
      CALL X04BAF(NERR,' ')
      WRITE (REC,FMT=99998) N
      DO 60 I = 1, 2
         CALL X04BAF(NERR,REC(I))
   60 CONTINUE
      WRITE (REC,FMT=99996) LWO, NMAXW, LIWORK, NMAXIW
      DO 80 I = 1, 2
         CALL X04BAF(NERR,REC(I))
   80 CONTINUE
      GO TO 220
  100 IF (MP.EQ.0) GO TO 160
      CALL X04BAF(NADV,' ')
      WRITE (REC,FMT=99997) NMAX, NMAXOK
      DO 120 I = 1, 3
         CALL X04BAF(NADV,REC(I))
  120 CONTINUE
      WRITE (REC,FMT=99996) LWO, NMAXW, LIWORK, NMAXIW
      DO 140 I = 1, 2
         CALL X04BAF(NADV,REC(I))
  140 CONTINUE
  160 MTNMAX = M*NMAXOK
      MTTWO = 2*M
      M4 = M + 1
      M5 = M4 + NMAXOK
      MT = MTNMAX*M
      M6 = M5 + MT
      M7 = M6 + MT
      M8 = M7 + MT
      M9 = M8 + MTNMAX
      M10 = M9 + MTNMAX
      M11 = M10 + MTNMAX
      M12 = M11 + NMAXOK
      M14 = M12
      IF (LIN.GT.1) M14 = M12 + MTNMAX
      M15 = M14 + MTNMAX
      M16 = M15 + MTTWO*M
      M22 = M16 + MTNMAX
      M23 = M22 + M
      IM2 = 1 + MTNMAX
      IM3 = IM2 + MTNMAX
      IM4 = IM3 + NMAXOK
      IM5 = IM4 + M
      IF (LIN.EQ.1 .OR. LIN.EQ.3) GO TO 180
      IM6 = IM5 + MTTWO + 1
      IM7 = IM6 + M + 1
      GO TO 200
  180 IM5 = 1
      IM6 = 1
      IM7 = 1
  200 CALL D02RAY(M,NMAXOK,N,NUMBEG,NUMMIX,MTNMAX,MTTWO,A,B,TOL,X,Y,IY,
     *            ABT,FCN,G,FCNEP,FCNA,FCNB,JACOBE,JACOBG,JACEPS,JACGEP,
     *            A1,B1,GAM,A2,B2,WORK(1),C1,D1,WORK(M4),WORK(M5)
     *            ,WORK(M6),WORK(M7),WORK(M8),WORK(M9),WORK(M10)
     *            ,WORK(M11),WORK(M12),WORK(M14),WORK(M15),WORK(M16)
     *            ,IWORK(1),IWORK(IM2),IWORK(IM3),IWORK(IM4)
     *            ,DELEPS,LP,MP,INIT,LIN,WORK(M22),WORK(M23),IWORK(IM5)
     *            ,2*M+1,IWORK(IM6),M+1,IWORK(IM7),IFLAG)
  220 RETURN
C
99999 FORMAT (' ONE OF THE FOLLOWING PARAMETER VALUES IS ILLEGAL',/'  ',
     *  '     N =',I8,/'      NP =',I8,/'     MNP =',I8,/'  NUMBEG =',
     *  I8,/'  NUMMIX =',I8)
99998 FORMAT (' A WORKSPACE ARRAY IS TOO SMALL - THE VALUE OF MNP USED',
     *  ' BY',/'  THE ROUTINE MUST BE .GE. 32 AND .GE. NP =',I6)
99997 FORMAT (' D02RAF WARNING',/'  THE VALUE OF MNP USED BY THE ROUTI',
     *  'NE HAS BEEN REDUCED',/'  FROM ',I6,' TO ',I6)
99996 FORMAT ('  LWORK  =',I6,'  THIS RESTRICTS THE VALUE OF MNP TO ',
     *  I6,/'  LIWORK =',I6,'  THIS RESTRICTS THE VALUE OF MNP TO ',I6)
      END
