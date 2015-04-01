      SUBROUTINE E04NCF(MM,N,NCLIN,LDA,LDR,A,BL,BU,CVEC,ISTATE,KX,X,R,B,
     *                  ITER,OBJ,CLAMDA,IW,LENIW,W,LENW,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1060 (JUL 1993).
C     MARK 17 REVISED. IER-1572 (JUN 1995).
C
C     ******************************************************************
C     E04NCF  solves problems of the form
C
C           Minimize               F(x)
C              x
C                                 (  x )
C           subject to    bl  .le.(    ).ge.  bu,
C                                 ( Ax )
C
C     where  '  denotes the transpose of a column vector,  x  denotes
C     the n-vector of parameters and  F(x) is one of the following
C     functions..
C
C     FP =  None                    (find a feasible point).
C     LP =  c'x
C     QP1=        1/2 x'Rx           R  n times n, symmetric pos. def.
C     QP2=  c'x + 1/2 x'Rx           .  .   ..        ..       ..  ..
C     QP3=        1/2 x'R'Rx         R  m times n, upper triangular.
C     QP4=  c'x + 1/2 x'R'Rx         .  .   ..  .   ..      ...
C     LS1=        1/2 (b - Rx)'(b - Rx)   R  m times n, rectangular.
C     LS2=  c'x + 1/2 (b - Rx)'(b - Rx)   .  .   ..  .     ...
C     LS3=        1/2 (b - Rx)'(b - Rx)  R  m times n, upper triangular.
C     LS4=  c'x + 1/2 (b - Rx)'(b - Rx)  .  .   ..  .   ..      ...
C
C     The matrix  R  is entered as the two-dimensional array  R  (of row
C     dimension  LDR).  If  LDR = 0,  R  is not accessed.
C
C     The vector  c  is entered in the one-dimensional array  CVEC.
C
C     NCLIN  is the number of general linear constraints (rows of  A).
C     (NCLIN may be zero.)
C
C     The first  N  elements of  BL  and   BU  are lower and upper
C     bounds on the variables.  The next  NCLIN  elements are
C     lower and upper bounds on the general linear constraints.
C
C     The matrix  A  of coefficients in the general linear constraints
C     is entered as the two-dimensional array  A  (of dimension
C     LDA  by  N).  If  NCLIN = 0,  A  is not accessed.
C
C     The vector  x  must contain an initial estimate of the solution,
C     and will contain the computed solution on output.
C
C
C     Complete documentation for  E04NCF  is contained in
C     Report SOL 86-1,
C     Users Guide for LSSOL (Version 1.0), by P.E. Gill,
C     S. J. Hammarling, W. Murray, M.A. Saunders and M.H. Wright,
C     Department of Operations Research, Stanford University, Stanford,
C     California 94305.
C
C     Systems Optimization Laboratory, Stanford University.
C
C     Version 1.0  Dated 24-Jan-1986.
C     Version 1.01 Dated 30-Jun-1986.   Level-2 BLAS added.
C     Version 1.02 Dated 13-May-1988.   Level-2 matrix routines added.
C     Version 1.03 Dated 19-Jun-1989.   Some obscure bugs fixed.
C     Version 1.04 Dated 26-Aug-1991.   NRANK bug fixed.
C     Version 1.05 Dated 20-Sep-1992.   Output modified.
C                        20-Oct-1992.   Summary file included.
C                        12-Jul-1994.   Hessian option added.
C                                       Debug printing eliminated.
C
C     Copyright  1983--1994  Stanford University.
C     This software is not in the public domain. Its use is governed by
C     a license agreement with Stanford University. It is illegal to
C     make copies except as authorized by the license agreement.
C
C     This material may be reproduced by or for the U.S. Government
C     pursuant to the copyright license under DAR clause 7-104.9(a)
C     (1979 Mar).
C
C     This material is based upon work partially supported by the
C     National Science Foundation under Grants MCS-7926009 and
C     ECS-8312142; the Department of Energy Contract AM03-76SF00326,
C     PA No. DE-AT03-76ER72018; the Army Research Office Contract DAA29-
C     84-K-0156; and the Office of Naval Research Grant N00014-75-C-0267
C     ******************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04NCF')
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, POINT3, POINT8
      PARAMETER         (ZERO=0.0D+0,POINT3=3.3D-1,POINT8=0.8D+0)
      DOUBLE PRECISION  POINT9, ONE, HUNDRD
      PARAMETER         (POINT9=0.9D+0,ONE=1.0D+0,HUNDRD=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJ
      INTEGER           IFAIL, ITER, LDA, LDR, LENIW, LENW, MM, N, NCLIN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  CLAMDA(N+NCLIN), CVEC(*), R(LDR,*), W(LENW),
     *                  X(N)
      INTEGER           ISTATE(N+NCLIN), IW(LENIW), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, DTMAX,
     *                  DTMIN, EPSPT3, EPSPT5, EPSPT8, EPSPT9, TOLACT,
     *                  TOLFEA, TOLRNK
      INTEGER           IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LDT, LDZY, LENNAM, LFORMH, LINES1,
     *                  LINES2, LPROB, MSGLS, NCOLT, NN, NNCLIN, NOUT,
     *                  NPROB
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPSVLS(MXPARM), WMACH(15)
      INTEGER           IPADLS(19), IPSVLS(MXPARM), LOCLS(LENLS)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMIN, CONDMX, CTX, EPSMCH, ERRMAX, FEAMAX,
     *                  FEAMIN, ROWNRM, RTEPS, SSQ1, SUMINF, XNORM
      INTEGER           I, IANRMJ, INFO, INFORM, J, JINF, JMAX, JSAVE,
     *                  LANORM, LAX, LCQ, LFEATL, LGQ, LITOTL, LJ,
     *                  LKACTV, LPX, LRES, LRES0, LRLAM, LT, LWRK,
     *                  LWTINF, LWTOTL, LZY, M, MAXACT, MAXNZ, MINACT,
     *                  MINFXD, MSGLVL, MXFREE, NACT1, NACTIV, NARTIF,
     *                  NCTOTL, NERR, NERROR, NFREE, NGQ, NRANK, NREJTD,
     *                  NRES, NRZ, NUMINF, NZ
      LOGICAL           COLD, FACTRZ, LINOBJ, NAMED, ROWERR, UNITQ,
     *                  VERTEX
      CHARACTER*2       PRBTYP
      CHARACTER*11      TITLE
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM)
      INTEGER           IPRMLS(MXPARM)
      CHARACTER*8       NAMES(1)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           F06KLF, P01ABF
      EXTERNAL          DNRM2, F06KLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04NBS, E04NBW, E04NCG, E04NCH,
     *                  E04NCM, E04NCS, E04NCU, E04NCW, E04NCX, E04NCY,
     *                  E04NCZ, F01QDF, F01QFF, F06DFF, F06FBF, F06FLF,
     *                  X02ZAZ, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (MSGLS,MSGLVL)
C     .. Save statement ..
      SAVE              /AX02ZA/, /DE04NC/, /EE04NC/
C     .. Data statements ..
      DATA              TITLE/' *** E04NCF'/
C     .. Executable Statements ..
C
C     Set the machine-dependent constants.
C
      CALL X02ZAZ
C
      EPSMCH = WMACH(3)
      RTEPS = WMACH(4)
      NOUT = WMACH(11)
      NERR = WMACH(12)
C
      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9
C
      NAMED = .FALSE.
C
      INFORM = 0
      ITER = 0
C
      CONDMX = MAX(ONE/EPSPT5,HUNDRD)
C
      NCTOTL = N + NCLIN
C
C     Set the default values of the parameters.
C
      CALL E04NCS(MM,N,NCLIN,TITLE)
C
C     Set all the parameters determined by the input code.
C
      IF (LPROB.EQ.1) THEN
         PRBTYP = 'FP'
         M = 0
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB.EQ.2) THEN
         PRBTYP = 'LP'
         M = 0
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB.EQ.3) THEN
         PRBTYP = 'QP'
         M = MM
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB.EQ.4) THEN
         PRBTYP = 'QP'
         M = MM
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB.EQ.5) THEN
         PRBTYP = 'QP'
         M = MM
         LINOBJ = .FALSE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB.EQ.6) THEN
         PRBTYP = 'QP'
         M = MM
         LINOBJ = .TRUE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB.EQ.7) THEN
         PRBTYP = 'LS'
         M = MM
         LINOBJ = .FALSE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB.EQ.8) THEN
         PRBTYP = 'LS'
         M = MM
         LINOBJ = .TRUE.
         FACTRZ = .TRUE.
      ELSE IF (LPROB.EQ.9) THEN
         PRBTYP = 'LS'
         M = MM
         LINOBJ = .FALSE.
         FACTRZ = .FALSE.
      ELSE IF (LPROB.EQ.10) THEN
         PRBTYP = 'LS'
         M = MM
         LINOBJ = .TRUE.
         FACTRZ = .FALSE.
      END IF
C
C     Assign the dimensions of arrays in the parameter list of E04NCZ.
C     Economies of storage are possible if the minimum number of active
C     constraints and the minimum number of fixed variables are known in
C     advance.  The expert user should alter MINACT and MINFXD
C     accordingly.
C     If a linear program is being solved and the matrix of general
C     constraints is fat,  i.e.,  NCLIN .LT. N,  a non-zero value is
C     known for MINFXD.  Note that in this case, VERTEX must be
C     set  .TRUE..
C
      MINACT = 0
      MINFXD = 0
C
      VERTEX = .FALSE.
      IF ((PRBTYP.EQ.'LP' .OR. PRBTYP.EQ.'FP') .AND. NCLIN.LT.N) THEN
         MINFXD = N - NCLIN - 1
         VERTEX = .TRUE.
      END IF
C
      MXFREE = N - MINFXD
      MAXACT = MAX(1,MIN(N,NCLIN))
      MAXNZ = N - (MINFXD+MINACT)
C
      IF (NCLIN.EQ.0) THEN
         LDZY = 1
         LDT = 1
         NCOLT = 1
         VERTEX = .FALSE.
      ELSE
         LDZY = MAX(1,MXFREE)
         LDT = MAX(MAXNZ,MAXACT)
         NCOLT = MXFREE
      END IF
C
C     Allocate certain arrays that are not done in E04NCM.
C
      LITOTL = 0
C
      LAX = 1
      LWTOTL = LAX + NCLIN - 1
C
C     Allocate remaining work arrays.
C
      CALL E04NCM(LPROB,N,NCLIN,LITOTL,LWTOTL)
C
      COLD = LCRASH .EQ. 0
C
C     Check input parameters and storage limits.
C
      LENNAM = 1
C
      CALL E04NBS(NERROR,MSGLVL,LCRASH,( .NOT. FACTRZ),LENIW,LENW,
     *            LITOTL,LWTOTL,N,NCLIN,ISTATE,KX,NAMED,NAMES,BIGBND,BL,
     *            BU,LPROB,LDA,LDR,MM,NERR,IFAIL)
C
      IF (NERROR.GT.0) THEN
         INFORM = 6
         GO TO 100
      END IF
C
      LKACTV = LOCLS(1)
C
      LANORM = LOCLS(2)
      LPX = LOCLS(4)
      LRES = LOCLS(5)
      LRES0 = LOCLS(6)
      LGQ = LOCLS(8)
      LCQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LZY = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK = LOCLS(14)
      LFEATL = LOCLS(15)
C
      IF (TOLFEA.GT.ZERO) CALL F06FBF(N+NCLIN,TOLFEA,W(LFEATL),1)
C
      IANRMJ = LANORM
      DO 20 J = 1, NCLIN
         W(IANRMJ) = DNRM2(N,A(J,1),LDA)
         IANRMJ = IANRMJ + 1
   20 CONTINUE
      IF (NCLIN.GT.0) CALL F06FLF(NCLIN,W(LANORM),1,ASIZE,AMIN)
C
      CALL F06FLF(NCTOTL,W(LFEATL),1,FEAMAX,FEAMIN)
      CALL DCOPY(NCTOTL,W(LFEATL),1,W(LWTINF),1)
      CALL DSCAL(NCTOTL,(ONE/FEAMIN),W(LWTINF),1)
C
      SSQ1 = ZERO
C
      IF (FACTRZ) THEN
C        ===============================================================
C        Factorize R using QR or Cholesky.  KX must be initialized.
C        ===============================================================
         DO 40 I = 1, N
            KX(I) = I
   40    CONTINUE
C
         IF (PRBTYP.EQ.'LP' .OR. PRBTYP.EQ.'FP') THEN
            NRANK = 0
         ELSE IF (PRBTYP.EQ.'QP') THEN
C           ------------------------------------------------------------
C           Compute the Cholesky factorization of R.  The Hessian is
C           M times M and resides in the upper left-hand corner of R.
C           ------------------------------------------------------------
            DO 60 J = M + 1, N
               CALL F06FBF(M,ZERO,R(1,J),1)
   60       CONTINUE
C
            CALL E04NCW(LDR,M,NRANK,TOLRNK,KX,R,MSGLVL,INFO)
C
            IF (NRANK.GT.0) CALL F06FBF(NRANK,ZERO,W(LRES0),1)
C
         ELSE IF (PRBTYP.EQ.'LS') THEN
C           ------------------------------------------------------------
C           Compute the orthogonal factorization PRQ = ( U ),  where P
C                                                      ( 0 )
C           is an orthogonal matrix and Q is a permutation matrix.
C           Overwrite R with the upper-triangle U.  The orthogonal
C           matrix P is applied to the residual and discarded.  The
C           permutation is stored in the array KX.  Once U has been
C           computed we need only work with vectors of length N within
C           E04NCZ.  However, it is necessary to store the sum of
C           squares of the terms  B(NRANK+1),...,B(M),  where B = Pr.
C           ------------------------------------------------------------
            CALL F01QFF('Column interchanges',M,N,R,LDR,W(LWRK),
     *                  IW(LKACTV),W(LGQ),INFO)
C
            LJ = LKACTV
            DO 80 J = 1, N
               JMAX = IW(LJ)
               IF (JMAX.GT.J) THEN
                  JSAVE = KX(JMAX)
                  KX(JMAX) = KX(J)
                  KX(J) = JSAVE
               END IF
               LJ = LJ + 1
   80       CONTINUE
C
            CALL F01QDF('Transpose','Separate',M,MIN(N,M-1),R,LDR,
     *                  W(LWRK),1,B,M,W(LGQ),INFO)
C
            ROWNRM = DNRM2(N,R(1,1),LDR)
            IF (ROWNRM.LE.TOLRNK .OR. ABS(R(1,1)).LE.ROWNRM*TOLRNK) THEN
               NRANK = 0
            ELSE
               NRANK = F06KLF(MIN(N,M),R,LDR+1,TOLRNK)
            END IF
C
            IF (M.GT.NRANK) SSQ1 = DNRM2(M-NRANK,B(NRANK+1),1)
C
            IF (NRANK.GT.0) CALL DCOPY(NRANK,B,1,W(LRES0),1)
         END IF
      ELSE
C        ===============================================================
C        R is input as an upper-triangular matrix with M rows.
C        ===============================================================
         NRANK = M
         IF (NRANK.GT.0) THEN
            IF (PRBTYP.EQ.'QP') THEN
               CALL F06FBF(NRANK,ZERO,W(LRES0),1)
            ELSE IF (PRBTYP.EQ.'LS') THEN
               CALL DCOPY(NRANK,B,1,W(LRES0),1)
            END IF
         END IF
      END IF
C
      IF (MSGLVL.GT.0 .AND. NRANK.LT.N .AND. PRBTYP.NE.'LP' .AND.
     *    PRBTYP.NE.'FP') THEN
         WRITE (REC,FMT=99988) NRANK
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
C     ------------------------------------------------------------------
C     Find an initial working set.
C     ------------------------------------------------------------------
      CALL E04NCU(COLD,VERTEX,NCLIN,NCTOTL,NACTIV,NARTIF,NFREE,N,LDA,
     *            ISTATE,IW(LKACTV),BIGBND,TOLACT,A,W(LAX),BL,BU,X,
     *            W(LGQ),W(LWRK))
C
C     ------------------------------------------------------------------
C     Compute the TQ factorization of the constraints while keeping R in
C     upper-triangular form.  Transformations associated with Q are
C     applied to CQ.  Transformations associated with P are applied to
C     RES0.  If some simple bounds are in the working set,  KX is
C     re-ordered so that the free variables come first.
C     ------------------------------------------------------------------
C     First, add the bounds. To save a bit of work, CQ is not loaded
C     until after KX has been re-ordered.
C
      NGQ = 0
      NRES = 0
      IF (NRANK.GT.0) NRES = 1
      UNITQ = .TRUE.
C
      CALL E04NCX(UNITQ,INFORM,NZ,NFREE,NRANK,NRES,NGQ,N,LDZY,LDA,LDR,
     *            LDT,ISTATE,KX,CONDMX,A,R,W(LT),W(LRES0),W(LCQ),W(LZY),
     *            W(LWRK),W(LPX),W(LRLAM),MSGLVL)
C
      IF (LINOBJ) THEN
C
C        Install the transformed linear term in CQ.
C        E04NBW applies the permutations in KX to CVEC.
C
         NGQ = 1
         CALL DCOPY(N,CVEC,1,W(LCQ),1)
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,W(LCQ),W(LZY),W(LWRK))
      END IF
C
      IF (NACTIV.GT.0) THEN
         NACT1 = NACTIV
         NACTIV = 0
C
         CALL E04NCY(UNITQ,VERTEX,INFORM,1,NACT1,NACTIV,NARTIF,NZ,NFREE,
     *               NRANK,NREJTD,NRES,NGQ,N,LDZY,LDA,LDR,LDT,ISTATE,
     *               IW(LKACTV),KX,CONDMX,A,R,W(LT),W(LRES0),W(LCQ),
     *               W(LZY),W(LWRK),W(LPX),W(LRLAM),MSGLVL)
      END IF
C
C     ------------------------------------------------------------------
C     Move the initial  x  onto the constraints in the working set.
C     Compute the transformed residual vector  Pr = Pb - RQ'x.
C     ------------------------------------------------------------------
      CALL E04NCH(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NRANK,NZ,N,
     *            NCTOTL,LDZY,LDA,LDR,LDT,ISTATE,IW(LKACTV),KX,JMAX,
     *            ERRMAX,CTX,XNORM,A,W(LAX),BL,BU,W(LCQ),W(LRES),
     *            W(LRES0),W(LFEATL),R,W(LT),X,W(LZY),W(LWRK),W(LPX))
C
      IF (ROWERR) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99981)
            CALL X04BAF(IPRINT,REC(1))
         END IF
         INFORM = 3
         NUMINF = 1
         SUMINF = ERRMAX
         GO TO 100
      END IF
C
      JINF = 0
C
      CALL E04NCZ(PRBTYP,NAMED,NAMES,LINOBJ,UNITQ,INFORM,ITER,JINF,
     *            NCLIN,NCTOTL,NACTIV,NFREE,NRANK,NZ,NRZ,N,LDA,LDR,
     *            ISTATE,IW(LKACTV),KX,CTX,OBJ,SSQ1,SUMINF,NUMINF,XNORM,
     *            BL,BU,A,CLAMDA,W(LAX),W(LFEATL),R,X,W)
C
C     ==================================================================
C     If required, form the triangular factor of the Hessian.
C     ==================================================================
      IF (LFORMH.GT.0) THEN
         CALL E04NCG('Permuted Hessian',UNITQ,NFREE,N,NRANK,LDZY,LDR,KX,
     *               R,W(LZY),W(LWRK),W(LPX))
      END IF
C
      OBJ = OBJ + CTX
      IF (PRBTYP.EQ.'LS' .AND. NRANK.GT.0) CALL DCOPY(NRANK,W(LRES),1,B,
     *    1)
C
C     ==================================================================
C     Print messages if required.
C     ==================================================================
  100 IF (MSGLVL.GT.0) THEN
         IF (INFORM.EQ.0) THEN
            IF (PRBTYP.EQ.'FP') THEN
               WRITE (REC,FMT=99999)
            ELSE
               WRITE (REC,FMT=99998) PRBTYP
            END IF
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (INFORM.EQ.1) WRITE (REC,FMT=99997) PRBTYP
         IF (INFORM.EQ.2) WRITE (REC,FMT=99996) PRBTYP
         IF (INFORM.EQ.3) WRITE (REC,FMT=99995)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99994)
         IF (INFORM.EQ.5) WRITE (REC,FMT=99993)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99992) NERROR
         IF (INFORM.GE.1 .AND. INFORM.LE.6) CALL X04BAY(IPRINT,2,REC)
C
         IF (INFORM.LT.6) THEN
            IF (NUMINF.EQ.0) THEN
               IF (PRBTYP.NE.'FP') THEN
                  WRITE (REC,FMT=99991) PRBTYP, OBJ
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            ELSE IF (INFORM.EQ.3) THEN
               WRITE (REC,FMT=99990) SUMINF
               CALL X04BAY(IPRINT,2,REC)
            ELSE
               WRITE (REC,FMT=99989) SUMINF
               CALL X04BAY(IPRINT,2,REC)
            END IF
            IF (NUMINF.GT.0) OBJ = SUMINF
         END IF
      END IF
C
C     Recover the optional parameters set by the user.
C
      CALL F06DFF(MXPARM,IPSVLS,1,IPRMLS,1)
      CALL DCOPY(MXPARM,RPSVLS,1,RPRMLS,1)
C
      IF (INFORM.GT.0 .AND. (IFAIL.EQ.0 .OR. IFAIL.EQ.-1)) THEN
         IF (INFORM.EQ.1) WRITE (REC,FMT=99987) PRBTYP
         IF (INFORM.EQ.2) WRITE (REC,FMT=99986) PRBTYP
         IF (INFORM.EQ.3) WRITE (REC,FMT=99985)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99984)
         IF (INFORM.EQ.5) WRITE (REC,FMT=99983)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99982) NERROR
         CALL X04BAY(NERR,2,REC)
      END IF
      IFAIL = P01ABF(IFAIL,INFORM,SRNAME,0,REC)
      RETURN
C
C
C
C     End of  E04NCF. (LSSOL)
C
99999 FORMAT (/' Exit E04NCF - Feasible point found.')
99998 FORMAT (/' Exit E04NCF - Optimal ',A2,' solution.')
99997 FORMAT (/' Exit E04NCF - Weak ',A2,' solution.')
99996 FORMAT (/' Exit E04NCF - ',A2,' solution is unbounded.')
99995 FORMAT (/' Exit E04NCF - Cannot satisfy the linear constraints. ')
99994 FORMAT (/' Exit E04NCF - Too many iterations.')
99993 FORMAT (/' Exit E04NCF - Too many iterations without changing X.')
99992 FORMAT (/' Exit E04NCF - ',I7,' errors found in the input parame',
     *       'ters.  Problem abandoned.')
99991 FORMAT (/' Final ',A2,' objective value =',G16.7)
99990 FORMAT (/' Minimum sum of infeasibilities =',G16.7)
99989 FORMAT (/' Final sum of infeasibilities =',G16.7)
99988 FORMAT (/' Rank of the objective function data matrix = ',I5)
99987 FORMAT (/' ** Weak ',A2,' solution.')
99986 FORMAT (/' ** ',A2,' solution is unbounded.')
99985 FORMAT (/' ** Cannot satisfy the linear constraints. ')
99984 FORMAT (/' ** Too many iterations.')
99983 FORMAT (/' ** Too many iterations without changing X.')
99982 FORMAT (/' ** ',I7,' errors found in the input parameters.  Prob',
     *       'lem abandoned.')
99981 FORMAT (' XXX  Cannot satisfy the working set constraints to the',
     *       ' accuracy requested.')
      END
