        SUBROUTINE CCHDEC(B,A,N,IER)
 
C       CHOLESKY DECOMPOSITION OF COMPLEX HERMITIAN MATRIX; B IS INPUT
C       MATRIX OF ORDER N IN SYMMETRIC STORAGE MODE (REQUIRES ARRAY OF
C       LENGTH N*(N+1)/2 IN CALLING ROUTINE); A IS OUTPUT (IN SAME FORM)
C       OF LOWER TRIANGULAR MATRIX
 
        COMPLEX A(*),B(*),T
 
        IER=0
 
        DO 30 I=1,N
        II=(I-1)*I/2
 
           DO 15 J=1,I-1
           JJ=(J-1)*J/2
           T=B(II+J)
           I1=II
 
              DO 10 K=1,J-1
              I1=I1+1
              JJ=JJ+1
10            T=T-A(I1)*CONJG(A(JJ))
 
15         A(II+J)=T/CONJG(A(JJ+1))
 
        I1=II
        T=B(II+I)
 
           DO 20 K=1,I-1
           I1=I1+1
20         T=T-A(I1)*CONJG(A(I1))
 
        IF (REAL(T).LT.0) THEN
           IER=-1
           RETURN
           END IF
 
30      A(II+I)=CMPLX(SQRT(REAL(T)),0.)
 
        RETURN
        END
C
C**************************
C
        SUBROUTINE CCHDECI(B,N,IER)
 
C       CHOLESKY DECOMPOSITION OF COMPLEX HERMITIAN MATRIX; B IS INPUT
C       MATRIX OF ORDER N IN SYMMETRIC STORAGE MODE (REQUIRES ARRAY OF
C       LENGTH N*(N+1)/2 IN CALLING ROUTINE);
c        in this version output overwrites input 

        COMPLEX B(*),T
 
        IER=0
 
        DO 30 I=1,N
        II=(I-1)*I/2
 
           DO 15 J=1,I-1
           JJ=(J-1)*J/2
           T=B(II+J)
           I1=II
 
              DO 10 K=1,J-1
              I1=I1+1
              JJ=JJ+1
10            T=T-B(I1)*CONJG(B(JJ))
 
15         B(II+J)=T/CONJG(B(JJ+1))
 
        I1=II
        T=B(II+I)
 
           DO 20 K=1,I-1
           I1=I1+1
20         T=T-B(I1)*CONJG(B(I1))
 
        IF (REAL(T).LT.0) THEN
           IER=-1
           RETURN
           END IF
 
30      B(II+I)=CMPLX(SQRT(REAL(T)),0.)
 
        RETURN
        END
C
C**************************
C
        SUBROUTINE CLTSLV(A,N,B,M)
 
C       LINEAR EQUATION SOLVER FOR COMPLEX LOWER TRIANGULAR SYSTEM
C       THE EQUATION TO SOLVE IS AX=B WHERE A IS N x N AND
C       B IS N x M; A IS LOWER TRIANGULAR IN SYMMETRIC STORAGE MODE
C       RESULT IS RETURNED IN B
 
        COMPLEX A(*),B(N,M),TEMP
 
        DO 50 I=1,M
 
           DO 20 J=1,N
 
           J1=J*(J-1)/2+1
           JL=J1+J-2
 
             TEMP=0.0
             DO 10 K=J1,JL
             TEMP=TEMP+A(K)*B(K-J1+1,I)
10           CONTINUE
 
           B(J,I)=(B(J,I)-TEMP)/A(JL+1)
20         CONTINUE
50      CONTINUE
        RETURN
        END
c__________________________________________________________________
c
      subroutine cutslv(a,n,b,m)

c       linear equation solver for complex upper triangular system
c       the equation to solve is atx=b where a is n x n and
c       b is n x m; a is lower triangular in symmetric storage mode
c       (and at is Hermitian transpose)
c       result is returned in b

      complex a(*),b(n,m)
      integer j,i,k,n,m,kk

      do j=1,n
         do k = 1,j-1
            kk = ((n-k)*(n-k+1))/2 + n-j+1
            do i = 1,m
                 b(n-j+1,i) = b(n-j+1,i)-conjg(a(kk))*b(n-k+1,i)
            enddo
         enddo
         kk = ((n-j+1)*(n-j+2))/2
         do i = 1,m
            b(n-j+1,i) = b(n-j+1,i)/conjg(a(kk))
         enddo
      enddo
      return
      end
c
c************************************************
c
       subroutine cinv2(a)
c       quick and dirty complex 2 x2 inverse
       complex a(2,2),det,temp
       ierr = 0
       det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
       temp = a(1,1)/det
       a(1,1) = a(2,2)/det
       a(2,2) = temp
       a(2,1) = -a(2,1)/det
       a(1,2) = -a(1,2)/det
       return
       end
C
C********************************************
C
        SUBROUTINE CMMULT(A,B,C,NN,NP,NQ)
 
C       MULTIPLIES N x P MATRIX A BY P x Q MATRIX B AND PUTS RESULT IN C
 
        COMPLEX A(NN,NP),B(NP,NQ),C(NN,NQ),TEMP
 
        DO 10 IN=1,NN
        DO 10 IQ=1,NQ
        TEMP= (0.,0.)
 
           DO  5 IP=1,NP
5          TEMP=TEMP+A(IN,IP)*B(IP,IQ)
 
10      C(IN,IQ)=TEMP
        RETURN
        END
c______________________________________________________________________
c
        subroutine cmus1(s,n,m,u,v)
c       pre-multiplies n x n complex hermitian matrix s (in symmetric storage
c        mode) by m x n complex matrix u; output is in v
 
        complex s(*), u(m,n),v(m,n),temp
 
        do 20 j = 1,n
        do 20 i = 1,m
        temp = (0.,0.)
           do 10 k = 1,n
           if(k.ge.j) then
              ii = (k*(k-1))/2+j
              temp = temp + u(i,k)*s(ii)
           else
              ii = (j*(j-1))/2+k
              temp = temp + u(i,k)*conjg(s(ii))
           end if
10         continue
        v(i,j) = temp
20      continue
        return
        end
c___________________________________________________________
c
      subroutine usuc1(s,u,n,nk,work,t)
c
c       complex multiplication routine; computes u s u^
c      result (output in t) may overwrite s

       complex s(*),t(*),u(nk,n),work(nk,n),temp
 
       call cmus1(s,n,nk,u,work)
 
       ii = 0
          do 20 i = 1,nk
          do 20 j = 1,i
             ii = ii + 1
             temp  = (0.,0.)
             do 10 k = 1,n
10           temp = temp + conjg(u(j,k))* work(i,k)
 
20       t(ii) = temp
      return
      end
c
c*******************************************
c
      subroutine mkarai(s,n,ar,ai)
      complex s(*)
      real ar(n,n),ai(n,n)

      ij = 0
         do 10 i = 1,n
            do 5 j = 1,i-1
            ij = ij + 1
            ar(i,j) = real(s(ij))
            ar(j,i) = real(s(ij))
            ai(i,j) = aimag(s(ij))
5           ai(j,i) = -aimag(s(ij))
         ij = ij + 1
         ar(i,i) = real(s(ij))
10       ai(i,i) = 0.0
       return
       end
c______________________________________________________________________
c
      SUBROUTINE HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
 
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL AR(NM,N),AI(NM,N),D(N),E(N),E2(N),TAU(2,N)
      REAL F,G,H,FI,GI,HH,SI,SCALE

      TAU(1,N) = 1.0
      TAU(2,N) = 0.0

      DO 100 I = 1,N
100     D(I) = AR(I,I)
********** FOR I=N STEP -1  UNTIL 1 DO -- *********
        DO 300 II=1,N
           I=N+1-II
           L=I-1
           H=0.0
           SCALE=0.0
           IF(L .LT. 1) GO TO 130
C********* SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
           DO 120 K=1,L
120        SCALE=SCALE+ABS(AR(I,K)) + ABS(AI(I,K))
 
           IF (SCALE .NE. 0.0) GO TO 140
           TAU(1,L) = 1.0
           TAU(2,L) = 0.0
130        E(I)=0.0
           E2(I)=0.0
           GO TO 290
 
140        DO 150 K=1,L
              AR(I,K)=AR(I,K)/SCALE
              AI(I,K)=AI(I,K)/SCALE
              H=H+AR(I,K)*AR(I,K)+AI(I,K)*AI(I,K)
150        CONTINUE
 
           E2(I)=SCALE * SCALE * H
           G = SQRT(H)
           E(I) = SCALE * G
           F = CABS(CMPLX(AR(I,L),AI(I,L)))
C************* FORM NEXT DIAGONAL EOLEMENT OF MATRIX T
           IF(F .EQ. 0.0) GO TO 160
           TAU(1,L)=(AI(I,L) * TAU(2,I)-AR(I,L)*TAU(1,I))/F
           SI=(AR(I,L)*TAU(2,I)+AI(I,L)*TAU(1,I))/F
           H=H+F*G
           G=1.0+G/F
           AR(I,L)=G*AR(I,L)
           AI(I,L)=G*AI(I,L)
           IF(L.EQ.1) GO TO 270
           GO TO 170
160        TAU(1,L)=-TAU(1,I)
           SI = TAU(2,I)
           AR(I,L) = G
170        F=0.0
 
           DO 240 J=1,L
              G=0.0
              GI=0.0
C***** FORM ELEMENT OF A*U  ***********
              DO 180 K=1,J
                 G=G+AR(J,K)*AR(I,K)+AI(J,K)*AI(I,K)
                 GI = GI - AR(J,K)*AI(I,K)+AI(J,K)*AR(I,K)
180           CONTINUE
 
              JP1=J+1
              IF(L .LT. JP1) GO TO 220
 
              DO 200 K=JP1,L
                 G=G+AR(K,J)*AR(I,K)-AI(K,J)*AI(I,K)
                 GI=GI-AR(K,J)*AI(I,K)-AI(K,J)*AR(I,K)
200           CONTINUE
C************FORM ELEMENT OF P*****************
220           E(J)=G / H
              TAU(2,J)=GI/H
              F=F+E(J)*AR(I,J)-TAU(2,J)*AI(I,J)
240        CONTINUE
 
           HH=F/(H+H)
C ************ FORM REDUCED A **************
           DO 260 J=1,L
              F=AR(I,J)
              G=E(J)-HH*F
              E(J)=G
              FI=-AI(I,J)
              GI=TAU(2,J)-HH*FI
              TAU(2,J) = -GI
 
              DO 260 K=1,J
                 AR(J,K)=AR(J,K)-F*E(K)-G*AR(I,K)+FI*TAU(2,K)+GI*
     1               AI(I,K)
                 AI(J,K)=AI(J,K)-F*TAU(2,K)-G*AI(I,K)-FI*E(K)-GI*
     1               AR(I,K)
260        CONTINUE
 
270        DO 280 K=1,L
              AR(I,K)=SCALE*AR(I,K)
              AI(I,K)=SCALE*AI(I,K)
280        CONTINUE
 
           TAU(2,L)=-SI
290        HH=D(I)
           D(I)=AR(I,I)
           AR(I,I)=HH
           AI(I,I)=SCALE*SQRT(H)
300     CONTINUE
 
        RETURN
        END
 
c>>>>>>>>>>>>>>>>The following text comes from EG_TQL2                          
        SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
 
                INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR
                REAL D(N),E(N),Z(NM,N)
                REAL B,C,F,G,H,P,R,S,MACHEP,X,Y
 
C ******* MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C       THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC
 
        MACHEP = 1.
        X = 1.
        DO 10 I = 1,100
        X = X/10.
        Y = 1. - X
        IF (Y.EQ.1) GO TO 20
        MACHEP = MACHEP/10.
10      CONTINUE
20      CONTINUE

ccc          print*,'MACHEP',MACHEP
 
C       MACHEP = 1.2E-7
 
        IERR =0
        IF (N .EQ. 1) GO TO 1001
 
        DO 100 I=2,N
100     E(I-1)=E(I)
 
        F=0.0
        B=0.0
        E(N)=0.0
 
        DO 240 L=1,N
           J=0
           H=MACHEP*(ABS(D(L))+ABS(E(L)))
           IF(B .LT. H) B=H
C ********* LOOK FOR SMALL SUBDIAGONAL ELEMENT **********
           DO 110 M=L,N
              IF(ABS(E(M)) .LE. B) GO TO 120
C ******* E(N) IS ALWAYS ZERO SO THERE IS NO EXIT
C               THROUGH THE BOTTOM OF THE LOOP
110        CONTINUE
 
120        IF(M .EQ. L) GO TO 220
130        IF(J.EQ.30) GO TO 1000
           J=J+1
C ************* FORM SHIFT ************
           L1=L+1
           G=D(L)
           P=(D(L1)-G)/(2.0*E(L))
           R=SQRT(P*P+1.0)
           D(L)=E(L)/(P+SIGN(R,P))
           H=G-D(L)
 
           DO 140 I=L1,N
140        D(I)=D(I)-H
 
           F=F+H
C************ QL TRANSFORMATION *********
           P=D(M)
           C=1.0
           S=0.0
           MML=M-L
C************** FOR I = M-1 STEP -1 UNTIL L DO -- ***********
           DO 200 II=1,MML
              I=M-II
              G=C*E(I)
              H=C*P
              IF(ABS(P).LT.ABS(E(I))) GO TO 150
              C=E(I)/P
              R=SQRT(C*C+1.0)
              E(I+1)=S*P*R
              S=C/R
              C=1.0/R
              GO TO 160
150           C=P/E(I)
              R=SQRT(C*C+1.0)
              E(I+1)=S*E(I)*R
              S=1.0/R
              C=C*S
160           P=C*D(I)-S*G
              D(I+1)=H+S*(C*G+S*D(I))
C ********* FORM VECTOR ******
              DO 180 K=1,N
                 H=Z(K,I+1)
                 Z(K,I+1) = S*Z(K,I)+C*H
                 Z(K,I)=C*Z(K,I)-S*H
180           CONTINUE
 
200        CONTINUE
 
            E(L)=S*P
            D(L)=C*P
            IF(ABS(E(L)).GT.B) GO TO 130
220         D(L)=D(L)+F
240     CONTINUE
C************* ORDER EIGENVALUES AND EIGENVECTORS ********
        DO 300 II=2,N
           I=II-1
           K=I
           P=D(I)
 
           DO 260 J=II,N
              IF(D(J) .GE. P) GO TO 260
              K=J
              P=D(J)
260        CONTINUE
 
        IF(K .EQ. I) GO TO 300
        D(K)=D(I)
        D(I)=P
 
        DO 280 J=1,N
           P=Z(J,I)
           Z(J,I)=Z(J,K)
           Z(J,K)=P
280     CONTINUE
 
300     CONTINUE
 
           GO TO 1001
C  ************** SET ERROR -- NO CONVERGENCE TO AN EIGENVALUE AFTER
C                     30 ITERATIONS *************
1000    IERR = L
1001    RETURN
        END
c>>>>>>>>>>>>>>>>The following text comes from EG_HTRIBK                        
        SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
 
        INTEGER I,J,K,L,M,N,NM
        REAL AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
        REAL H,S,SI
 
        IF (M.EQ.0) GO TO 200
C********** TRANSFORM, EIGENVECTORS OF THE REAL SYMMETRIC
C          TRIDIAGANOL MTRIX TO THOSE OF THE HERMITIAN TRIDIAGONAL
C          MATRIX
        DO 50 K=1,N
 
            DO 50 J=1,M
               ZI(K,J)=-ZR(K,J)*TAU(2,K)
               ZR(K,J)=ZR(K,J)*TAU(1,K)
50      CONTINUE
 
        IF (N.EQ.1) GO TO 200
C********* RECOVER AND APPLY THE HOUSEHOLDER MATRICES
 
        DO 140 I=2,N
           L=I-1
           H=AI(I,I)
           IF(H.EQ. 0.0) GO TO 140
 
           DO 130 J=1,M
              S=0.0
              SI=0.0
 
              DO 110 K=1,L
                 S = S + AR(I,K) * ZR(K,J) -  AI(I,K) * ZI(K,J)
                 SI =  SI + AR(I,K) * ZI (K,J) + AI(I,K) * ZR(K,J)
110           CONTINUE
C********** Double divisions avoid posssible underflow******
              S = (S/H) / H
              SI = (SI / H)/ H
 
              DO 120 K = 1,L
                 ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
                 ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
120           CONTINUE
 
130        CONTINUE
 
140     CONTINUE
 
200     RETURN
        END
 
C
C*********************************
C
        SUBROUTINE CHEIG(N,AR,AI,FV1,FM1,ZR,ZI,W,IERR)
 
C       Computes eigenvalues and eigenvectors of Complex Hermitian
C       matrix
c
c       input: AR(N,N),AI(N,N) real and imaginary parts of input matrix
c              of order N
c
c              FV1(N),FM1(2,N) work space
c
c       output: W eigenvalues in ascending order
c
c              ZR,ZI real and imaginary parts of corresponding
c               (orthonormal) eigenvectors
c
c       Subroutine calls EISPACK routines HTRIDI,TQL2, and HTRIBK
c
 
        REAL AR(N,N),AI(N,N),ZR(N,N),ZI(N,N),W(N),FV1(N),FM1(2,N)
 
        CALL HTRIDI(N,N,AR,AI,W,FV1,ZR,FM1)
 
        DO 100 I=1,N
           DO 50 J=1,N
           ZR(I,J)=0.0
50         CONTINUE
        ZR(I,I)=1.0
100     CONTINUE

        CALL TQL2(N,N,W,FV1,ZR,IERR)
 
        IF (IERR .NE. 0) GO TO 99999
 
        CALL HTRIBK(N,N,AR,AI,FM1,N,ZR,ZI)
 
        GO TO 200
 
99999   WRITE(6,*) 'ERRROR IN C_H_EIG; ERROR NO. : ',IERR
 
200     RETURN
        END
c________________________________________________________________
c
      function cdot(u,v,n)
c         complex inner product
      complex u(n),v(n),cdot
      cdot = (0.,0.)
         do 10 i = 1,n
10       cdot = cdot + conjg(u(i))*v(i)
      return
      end 
c______________________________________________________________________
c
      subroutine sqrdc(x,ldx,n,p,qraux,jpvt,work,job)
      integer ldx,n,p,job
      integer jpvt(1)
      real x(ldx,1),qraux(1),work(1)
c
c     sqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     on entry
c
c        x       real(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    real(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   real(p).
c                qraux contains further information required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     sqrdc uses the following functions and subprograms.
c
c     blas saxpy,sdot,sscal,sswap,snrm2
c     fortran abs,amax1,min0,sqrt
c
c     internal variables
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      real maxnrm,snrm2,tt
      real sdot,nrmxl,t
      logical negj,swapj
c
c
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call sswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call sswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = snrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0e0
            maxj = l
            do 100 j = l, pu
               if (qraux(j) .le. maxnrm) go to 90
                  maxnrm = qraux(j)
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call sswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = 0.0e0
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = snrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0e0) go to 180
               if (x(l,l) .ne. 0.0e0) nrmxl = sign(nrmxl,x(l,l))
               call sscal(n-l+1,1.0e0/nrmxl,x(l,l),1)
               x(l,l) = 1.0e0 + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -sdot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call saxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (qraux(j) .eq. 0.0e0) go to 150
                     tt = 1.0e0 - (abs(x(l,j))/qraux(j))**2
                     tt = amax1(tt,0.0e0)
                     t = tt
                     tt = 1.0e0 + 0.05e0*tt*(qraux(j)/work(j))**2
                     if (tt .eq. 1.0e0) go to 130
                        qraux(j) = qraux(j)*sqrt(t)
                     go to 140
  130                continue
                        qraux(j) = snrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end
 
c>>>>>>>>>>>>>>>>the following text comes from eg_saxpy                         
      subroutine saxpy(n,sa,sx,incx,sy,incy)
c                                                                        sa20310
c     overwrite single precision sy with single precision sa*sx +sy.     sa20320
c     for i = 0 to n-1, replace  sy(ly+i*incy) with sa*sx(lx+i*incx) +
c       sy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n,
c       and ly is defined in a similar way using incy.
c                                                                        sa20330
      real sx(1),sy(1),sa                                                sa20340
      if(n.le.0.or.sa.eq.0.e0) return                                    sa20350
      if(incx.eq.incy) if(incx-1) 5,20,60                                sa20360
    5 continue                                                           sa20370
c                                                                        sa20380
c        code for nonequal or nonpositive increments.                    sa20390
c                                                                        sa20400
      ix = 1                                                             sa20410
      iy = 1                                                             sa20420
      if(incx.lt.0)ix = (-n+1)*incx + 1                                  sa20430
      if(incy.lt.0)iy = (-n+1)*incy + 1                                  sa20440
      do 10 i = 1,n                                                      sa20450
        sy(iy) = sy(iy) + sa*sx(ix)                                      sa20460
        ix = ix + incx                                                   sa20470
        iy = iy + incy                                                   sa20480
   10 continue                                                           sa20490
      return                                                             sa20500
c                                                                        sa20510
c        code for both increments equal to 1                             sa20520
c                                                                        sa20530
c                                                                        sa20540
c        clean-up loop so remaining vector length is a multiple of 4.    sa20550
c                                                                        sa20560
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40                                            sa20580
      do 30 i = 1,m                                                      sa20590
        sy(i) = sy(i) + sa*sx(i)                                         sa20600
   30 continue                                                           sa20610
      if( n .lt. 4 ) return                                              sa20620
   40 mp1 = m + 1                                                        sa20630
      do 50 i = mp1,n,4                                                  sa20640
        sy(i) = sy(i) + sa*sx(i)                                         sa20650
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)                             sa20660
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)                             sa20670
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)                             sa20680
   50 continue                                                           sa20690
      return                                                             sa20700
c                                                                        sa20710
c        code for equal, positive, nonunit increments.                   sa20720
c                                                                        sa20730
   60 continue                                                           sa20740
      ns = n*incx                                                        sa20750
          do 70 i=1,ns,incx                                              sa20760
          sy(i) = sa*sx(i) + sy(i)                                       sa20770
   70     continue                                                       sa20780
      return                                                             sa20790
      end                                                                sa20800
c>>>>>>>>>>>>>>>>the following text comes from eg_sdot                          
      real function sdot(n,sx,incx,sy,incy)
c                                                                        sd17440
c     returns the dot product of single precision sx and sy.             sd17450
c     sdot = sum for i = 0 to n-1 of  sx(lx+i*incx) * sy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c                                                                        sd17460
      real sx(1),sy(1)                                                   sd17470
      sdot = 0.0e0                                                       sd17480
      if(n.le.0)return                                                   sd17490
      if(incx.eq.incy) if(incx-1)5,20,60                                 sd17500
    5 continue                                                           sd17510
c                                                                        sd17520
c        code for unequal increments or nonpositive increments.          sd17530
c                                                                        sd17540
      ix = 1                                                             sd17550
      iy = 1                                                             sd17560
      if(incx.lt.0)ix = (-n+1)*incx + 1                                  sd17570
      if(incy.lt.0)iy = (-n+1)*incy + 1                                  sd17580
      do 10 i = 1,n                                                      sd17590
        sdot = sdot + sx(ix)*sy(iy)                                      sd17600
        ix = ix + incx                                                   sd17610
        iy = iy + incy                                                   sd17620
   10 continue                                                           sd17630
      return                                                             sd17640
c                                                                        sd17650
c        code for both increments equal to 1                             sd17660
c                                                                        sd17670
c                                                                        sd17680
c        clean-up loop so remaining vector length is a multiple of 5.    sd17690
c                                                                        sd17700
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40                                            sd17720
      do 30 i = 1,m                                                      sd17730
        sdot = sdot + sx(i)*sy(i)                                        sd17740
   30 continue                                                           sd17750
      if( n .lt. 5 ) return                                              sd17760
   40 mp1 = m + 1                                                        sd17770
      do 50 i = mp1,n,5                                                  sd17780
        sdot = sdot + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +                sd17790
     1   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4) sd17800
   50 continue                                                           sd17810
      return                                                             sd17820
c                                                                        sd17830
c        code for positive equal increments .ne.1.                       sd17840
c                                                                        sd17850
   60 continue                                                           sd17860
      ns=n*incx                                                          sd17870
      do 70 i=1,ns,incx                                                  sd17880
        sdot = sdot + sx(i)*sy(i)                                        sd17890
   70   continue                                                         sd17900
      return                                                             sd17910
      end                                                                sd17920
c>>>>>>>>>>>>>>>>the following text comes from eg_sscal                         
      subroutine sscal(n,sa,sx,incx)
c                                                                        ss36670
c     replace single precision sx by single precision sa*sx.             ss36680
c     for i = 0 to n-1, replace sx(1+i*incx) with  sa * sx(1+i*incx)
c                                                                        ss36690
      real sa,sx(1)                                                      ss36700
      if(n.le.0)return                                                   ss36710
      if(incx.eq.1)goto 20                                               ss36720
c                                                                        ss36730
c        code for increments not equal to 1.                             ss36740
c                                                                        ss36750
      ns = n*incx                                                        ss36760
          do 10 i = 1,ns,incx                                            ss36770
          sx(i) = sa*sx(i)                                               ss36780
   10     continue                                                       ss36790
      return                                                             ss36800
c                                                                        ss36810
c        code for increments equal to 1.                                 ss36820
c                                                                        ss36830
c                                                                        ss36840
c        clean-up loop so remaining vector length is a multiple of 5.    ss36850
c                                                                        ss36860
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40                                            ss36880
      do 30 i = 1,m                                                      ss36890
        sx(i) = sa*sx(i)                                                 ss36900
   30 continue                                                           ss36910
      if( n .lt. 5 ) return                                              ss36920
   40 mp1 = m + 1                                                        ss36930
      do 50 i = mp1,n,5                                                  ss36940
        sx(i) = sa*sx(i)                                                 ss36950
        sx(i + 1) = sa*sx(i + 1)                                         ss36960
        sx(i + 2) = sa*sx(i + 2)                                         ss36970
        sx(i + 3) = sa*sx(i + 3)                                         ss36980
        sx(i + 4) = sa*sx(i + 4)                                         ss36990
   50 continue                                                           ss37000
      return                                                             ss37010
      end                                                                ss37020
c>>>>>>>>>>>>>>>>the following text comes from eg_sswap                         
      subroutine sswap (n,sx,incx,sy,incy)
c                                                                        ss31870
c     interchange single precision sx and single precision sy.           ss31880
c     for i = 0 to n-1, interchange  sx(lx+i*incx) and sy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c                                                                        ss31890
      real sx(1),sy(1),stemp1,stemp2,stemp3                              ss31900
      if(n.le.0)return                                                   ss31910
      if(incx.eq.incy) if(incx-1) 5,20,60                                ss31920
    5 continue                                                           ss31930
c                                                                        ss31940
c       code for unequal or nonpositive increments.                      ss31950
c                                                                        ss31960
      ix = 1                                                             ss31970
      iy = 1                                                             ss31980
      if(incx.lt.0)ix = (-n+1)*incx + 1                                  ss31990
      if(incy.lt.0)iy = (-n+1)*incy + 1                                  ss32000
      do 10 i = 1,n                                                      ss32010
        stemp1 = sx(ix)                                                  ss32020
        sx(ix) = sy(iy)                                                  ss32030
        sy(iy) = stemp1                                                  ss32040
        ix = ix + incx                                                   ss32050
        iy = iy + incy                                                   ss32060
   10 continue                                                           ss32070
      return                                                             ss32080
c                                                                        ss32090
c       code for both increments equal to 1                              ss32100
c                                                                        ss32110
c                                                                        ss32120
c       clean-up loop so remaining vector length is a multiple of 3.     ss32130
c                                                                        ss32140
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40                                            ss32160
      do 30 i = 1,m                                                      ss32170
        stemp1 = sx(i)                                                   ss32180
        sx(i) = sy(i)                                                    ss32190
        sy(i) = stemp1                                                   ss32200
   30 continue                                                           ss32210
      if( n .lt. 3 ) return                                              ss32220
   40 mp1 = m + 1                                                        ss32230
      do 50 i = mp1,n,3                                                  ss32240
        stemp1 = sx(i)                                                   ss32250
        stemp2 = sx(i+1)                                                 ss32260
        stemp3 = sx(i+2)                                                 ss32270
        sx(i) = sy(i)                                                    ss32280
        sx(i+1) = sy(i+1)                                                ss32290
        sx(i+2) = sy(i+2)                                                ss32300
        sy(i) = stemp1                                                   ss32310
        sy(i+1) = stemp2                                                 ss32320
        sy(i+2) = stemp3                                                 ss32330
   50 continue                                                           ss32340
      return                                                             ss32350
   60 continue                                                           ss32360
c                                                                        ss32370
c     code for equal, positive, nonunit increments.                      ss32380
c                                                                        ss32390
      ns = n*incx                                                        ss32400
        do 70 i=1,ns,incx                                                ss32410
        stemp1 = sx(i)                                                   ss32420
        sx(i) = sy(i)                                                    ss32430
        sy(i) = stemp1                                                   ss32440
   70   continue                                                         ss32450
      return                                                             ss32460
      end                                                                ss32470
c>>>>>>>>>>>>>>>>the following text comes from eg_snrm2                         
      real function snrm2 ( n, sx, incx)
      integer          next
      real   sx(1),  cutlo, cuthi, hitest, sum, xmax, zero, one
      data   zero, one /0.0e0, 1.0e0/
c
c     euclidean norm of the n-vector stored in sx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 4.441e-16,  1.304e19 /
c
      if(n .gt. 0) go to 10
         snrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( abs(sx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( sx(i) .eq. zero) go to 200
      if( abs(sx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( abs(sx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200
c
  115 sum = sum + (sx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(abs(sx(j)) .ge. hitest) go to 100
   95    sum = sum + sx(j)**2
      snrm2 = sqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      snrm2 = xmax * sqrt(sum)
  300 continue
      return
      end
c>>>>>>>>>>>>>>>>the following text comes from eg_scopy                         
      subroutine scopy(n,sx,incx,sy,incy)
c                                                                        sc30540
c     copy single precision sx to single precision sy.                   sc30550
c     for i = 0 to n-1, copy  sx(lx+i*incx) to sy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c                                                                        sc30560
      real sx(1),sy(1)                                                   sc30570
      if(n.le.0)return                                                   sc30580
      if(incx.eq.incy) if(incx-1) 5,20,60                                sc30590
    5 continue                                                           sc30600
c                                                                        sc30610
c        code for unequal or nonpositive increments.                     sc30620
c                                                                        sc30630
      ix = 1                                                             sc30640
      iy = 1                                                             sc30650
      if(incx.lt.0)ix = (-n+1)*incx + 1                                  sc30660
      if(incy.lt.0)iy = (-n+1)*incy + 1                                  sc30670
      do 10 i = 1,n                                                      sc30680
        sy(iy) = sx(ix)                                                  sc30690
        ix = ix + incx                                                   sc30700
        iy = iy + incy                                                   sc30710
   10 continue                                                           sc30720
      return                                                             sc30730
c                                                                        sc30740
c        code for both increments equal to 1                             sc30750
c                                                                        sc30760
c                                                                        sc30770
c        clean-up loop so remaining vector length is a multiple of 7.    sc30780
c                                                                        sc30790
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40                                            sc30810
      do 30 i = 1,m                                                      sc30820
        sy(i) = sx(i)                                                    sc30830
   30 continue                                                           sc30840
      if( n .lt. 7 ) return                                              sc30850
   40 mp1 = m + 1                                                        sc30860
      do 50 i = mp1,n,7                                                  sc30870
        sy(i) = sx(i)                                                    sc30880
        sy(i + 1) = sx(i + 1)                                            sc30890
        sy(i + 2) = sx(i + 2)                                            sc30900
        sy(i + 3) = sx(i + 3)                                            sc30910
        sy(i + 4) = sx(i + 4)                                            sc30920
        sy(i + 5) = sx(i + 5)                                            sc30930
        sy(i + 6) = sx(i + 6)                                            sc30940
   50 continue                                                           sc30950
      return                                                             sc30960
c                                                                        sc30970
c        code for equal, positive, nonunit increments.                   sc30980
c                                                                        sc30990
   60 continue                                                           sc31000
      ns = n*incx                                                        sc31010
          do 70 i=1,ns,incx                                              sc31020
          sy(i) = sx(i)                                                  sc31030
   70     continue                                                       sc31040
      return                                                             sc31050
      end                                                                sc31060
c>>>>>>>>>>>>>>>>the following text comes from eg_sqrsl                         
      subroutine sqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      integer ldx,n,k,job,info
      real x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),xb(1)
c
c     sqrsl applies the output of sqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to sqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  sqrdc produces a factored orthogonal matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      real(ldx,p).
c               x contains the output of sqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in sqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to sqrdc.
c
c        qraux  real(p).
c               qraux contains the auxiliary output from sqrdc.
c
c        y      real(n)
c               y contains an n-vector that is to be manipulated
c               by sqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.

c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     real(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    real(n).
c               qty contains trans(q)*y, if its computation has
c               been requested.  here trans(q) is the
c               transpose of the matrix q.
c
c        b      real(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in sqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into sqrdc.)
c
c        rsd    real(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     real(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call sqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     sqrsl uses the following functions and subprograms.
c
c     blas saxpy,scopy,sdot
c     fortran abs,min0,mod
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      real sdot,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0e0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = 0.0e0
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call scopy(n,y,1,qy,1)
         if (cqty) call scopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0e0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -sdot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call saxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute trans(q)*y.
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0e0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -sdot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call saxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call scopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call scopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call scopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0e0
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0e0
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0e0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call saxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0e0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -sdot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call saxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -sdot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call saxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end
c
c************************************
c
        subroutine lsv(x,y,n,ip,b,var,wrk)
 
c     this version returns variances of estimates (but not covariances)
        real x(n,ip),y(n),b(ip),var(ip),wrk(*)
 
c       compute q-r decomp.
 
        job = 0
 
        call sqrdc(x,n,n,ip,wrk,jpvt,work,job)
 
c       compute least squares solution
 
        job = 110
        call sqrsl(x,n,n,ip,wrk,y,dum,y,b,y,dum,job,info)
 
c   compute diagonal part of error covariance matrix

        sig = 0.
        do 10 i = 1,n
        sig = sig+y(i)*y(i)
10      continue
 
        m = n
        if(n.gt.ip) m = m-ip
        sig = sig/m

          do 50 i = 1,ip

             do 40 j = i,ip
             if(j.eq.i) then
                var(j) = 1.
             else
                var(j) = 0.
             end if

                do 20 k = i,j-1
                var(j) = var(j)-var(k)*x(k,j)
20              continue

             var(j) = var(j)/x(j,j)
40           continue

             var(i) = var(i)*var(i)
             do 50 j = i+1,ip
             var(i) = var(i) + var(j)*var(j)
50        continue

        do 60 i = 1,ip
60      var(i) = var(i)*sig

        return
        end
