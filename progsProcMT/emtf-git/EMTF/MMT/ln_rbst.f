c______________________________________________________________________
c
      subroutine mks1(s,i1,n1,s11)
      complex s(*),s11(*)
      integer sindex,i1(n1)

ccc     given index sets i1,i2 and a SDM S (symetric storage mode)
ccc     reorder components and output partitioned SDM  ::
ccc     After reordering
ccc                     S = | S11  S12 |
ccc                         | S21  S22 |
ccc     where S11 is SDM for components in index set I1, etc.
ccc        S11 (which is Hermitian, is returned in
ccc        symetric storage mode (only lower triangular half of matrix)
ccc      VARIANT ON MKS1S2, WHICH ONLY RETURNS S11

ccc    S11:

      ii = 0
      do i = 1,n1
         do j = 1,i
            ii = ii + 1
            s11(ii) = s(sindex(i1(i),i1(j)))
         enddo
      enddo
      return
      end

c______________________________________________________________________
c
      subroutine ln_rbst(xr,xc,xp,var,s,nt,nd,nsta,ih,sn,tf,tf_se,
     &    nch_grp,itmax_stk,itmax_ln)

ccc   robust version of LNOISE
ccc   routine to estimate local noise variance, using robust regression
ccc   of each channel on PCs of other channels
ccc    ==> XR is input data (not modified)
ccc    ==> XC is input "cleaned" data (not modified)
ccc    ==> XP is ouptut "predicted" data
ccc    ==> SN is array of estimated "local noise" covariance matrices

      include 'iosize.inc'
      parameter (npcmin = 2,npcmax = ntmx-npcmin,n2max=nchmx)
      parameter (nwork=npcmax*(ntfmax+2*n2max+npcmax+100+1))
      complex xr(nd,nt),xc(nd,nt),xp(nd,nt),tf(nt,nt)
      real var(nt),eval(ntmx),pctol
      complex s(*),sn(*),s11(nsmx),u(ntmx,ntmx),work(nwork)
      real tf_se(nt,nt)
      integer ih(*),i1(ntmx),i2(ntmx),sindex,nch_grp(*),ipt_sn,
     &    itmax_stk,n_it_sdm

ccc   first see if variance estimates have been provided
      if(var(1) .eq. 0.0 ) then
ccc      variances haven't been estimated yet; just compute l1 norm of each channel
         varsum = 0.0
         do j = 1,nt
            var(j) = 0.0
            do i = 1,nd
               var(j) = var(j) + abs(xc(i,j))
            enddo
            var(j) = (var(j)/nd)**2.
            varsum = varsum + var(j)
         enddo
ccc      also, set tol (criteria for choosing # of PCs to use to .001
         pctol = .001
ccc      and do robust stacking of SDM for each group ...        2-28-97
         n_it_sdm = itmax_stk
      else
ccc      have estimates of noise variances, so don't use noise eigenvectors 
ccc      for PC regression
         pctol = 1.
ccc      and don't do robust stacking of SDM  ...         2-28-97
         n_it_sdm = 0
      endif

ccc   normalize "cleaned" data before finding PCs
      varsum = varsum/(nt*1000.)
      do j = 1,nt
         scale = sqrt(var(j)+varsum)
c         write(*,*) 'j,var(j) [ln_rbst]',j,var(j)
         do i = 1,nd
            xc(i,j) = xc(i,j)/scale
         enddo
      enddo

ccc   make SDM
      call rbstk2(xc,nt,nd,s,xp,n_it_sdm)

ccc    ::outer loop:: estimate local noise covariance for each "station"
ccc          (more generally, for each group of channels)

      do ista = 1,nsta
         isn = ipt_sn(nch_grp,ista)

ccc     compute index set I1
         n1 = 0
         n2 = 0
         do i =  1,nt
            if( (i.lt.ih(ista)).or.(i.ge.ih(ista+1)) ) then
               n1 = n1 + 1
               i1(n1) = i
            else
               n2 = n2 + 1
               i2(n2) = i
            end if 
         enddo
         if(n2.eq.0) goto 100
ccc      no data channels in this grouping ... go to end of loop

ccc      else reorder cross-products and partition matrix
         call mks1(s,i1,n1,s11)
ccc      compute principal components of normalized S11
         call rsp(s11,n1,n1,u,eval,0,work,dum)
ccc      decide on number of PCs to keep

         npc = 0
c         write(*,*) 'evals',eval
         do i = 1,n1
            if(eval(i).ge.pctol) then
               npc = npc + 1
            endif
         enddo
         npc = max(npc,npcmin)
         npc = min(npc,npcmax)
         npc = min(npc,nd-1)
c         write(*,*) 'n1 = ',n1
c         write(*,*) 'n2 = ',n2
c         write(*,'(6hUsing ,i2,22h PCs for Prediction TF)') npc

ccc      form NPC linear combinations of data
         call lc_dat(xc,nd,nt,u,n1,npc,i1,work)
ccc      do robust regression on PCs
ccc      the integers iw# serve to partition the array WORK
ccc      into the different pieces used inside RBSTREG
ccc      WORK should be dimensioned to allow for NPC*(ND+2*N2+NPC+NB+1)
ccc       complex numbers;  NB is block size inside RBSTREG, currently
ccc       set to 100 (bigger than needed, but who cares? ... this is
ccc       tiny compared to ND)
         iw2 = npc*nd+1       ! start of array b
         iw3 = iw2+npc*n2     ! start of array bsv
         iw4 = iw3+npc*n2     ! start of array xxinv
         iw5 = iw4+npc*npc    ! start of array work (inside rbstreg)
c         write(18,*) 'STATION ',ista
         call rbstreg(work,xr(1,ih(ista)),xp(1,ih(ista)),nd,npc,n2,
     &      work(iw2),work(iw3),work(iw5),sn(isn)
     &     ,work(iw4),itmax_ln,icvg,0)

c         write(*,*) 'nsn,ista,sn',nsn,ista,(sn(k),k=isn,isn+20)

c         if(icvg.lt.0) print*,'did not converge'

         call atf(tf,tf_se,nt,n1,i1,n2,i2,npc,work(iw2),u,work(iw4)
     &       ,sn(isn),var,varsum)

c       write(*,*) 'after atf:nsn,ista,sn',nsn,ista,(sn(k),k=isn,isn+20)

100      continue
      enddo    ! do ista = 1,nsta
      return
      end
c______________________________________________________________________
c
      subroutine wrt_sn(iounit,sn,nch,nsta,nch_grp)
      integer nch(nsta),nsn,ipt_sn
      complex sn(*)
      include 'nstamx.inc'
      real sdiag(ntmx)

      write(iounit,*) 'LOCAL NOISE (Variance; Correlation'
      do ista = 1,nsta
         isn = ipt_sn(nch,ista)
c          extract diagonal (variances)
         call diag(sn(isn),nch(ista),sdiag)

         do k=1,nch(ista)
            sdiag(k) = sqrt(sdiag(k))
         enddo

         k1 = 0
         do i = 1,nch(ista)
          if(i.eq.1) then
             write(iounit,21) real(sn(isn)),1.0,0.0,ista
          else
             write(iounit,20) Real(sn(isn+k1+i-1)),
     &        ( sn(isn+k1+j-1)/(sdiag(i)*sdiag(j)),j=1,i)
           endif
           k1 = k1 + i
         enddo
      enddo
      return

c10    format(1x,5(e12.4,2x))
20    format(e12.4,2x,10(f5.3,f5.3,2x))
21    format(e12.4,2x,f5.3,f5.3,20x,i2)
      end
c______________________________________________________________________
c
      subroutine wrt_cormat(iounit,s,nt)
      complex s(*)
      real sdiag(100)
      integer iounit

      call diag(s,nt,sdiag)
c          write out variances
c      write(iounit,*) 'VARIANCES'
c      write(iounit,10) (sdiag(k),k=1,nt)
      do k=1,nt
         sdiag(k) = sqrt(sdiag(k))
      enddo

c          write out correlation matrix
c      write(iounit,*) 'Squared Correlations'
      k1 = 0
      do i = 1,nt
         write(iounit,30) ( abs(s(k1+j)/(sdiag(i)*sdiag(j)))**2,j=1,i)
         k1 = k1 + i
         if(mod(i,5).eq.0 ) write(iounit,*)
      enddo
      return

c10    format(5(e12.4,2x))
20    format(15(2f4.2,1x))
30    format(3(5f6.3,3x))
      end
c______________________________________________________________________
c
      subroutine lc_dat(x,nd,nt,u,n1,nlc,i1,pcx)
      complex u(n1,nlc),x(nd,nt),pcx(nd,nlc)
      integer i1(n1)

ccc      forms NLC linear combinations of N1 data channels
ccc       (which are part of a potentially larger number NT of
ccc       complex data channels).  Input data matrix is X;
ccc       ouput is in PCX

      do i = 1,nd
         do j = 1,nlc
            pcx(i,j) = (0.,0.)
            do k = 1,n1
               pcx(i,j) = pcx(i,j) + conjg(u(k,j))*x(i,i1(k))
            enddo
         enddo
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine atf(tf,tf_se,nt,n1,i1,n2,i2,npc,b,u,xxinv,sn,var,
     &   varsum)
ccc      copy TFs into array tf
ccc      NOTE: VAR is variances used to scale cleaned data
ccc          before constructing PCs; SN is noise covariance
ccc          estimated by rbstreg
      complex tf(nt,nt),b(npc,n2),u(n1,npc),xxinv(npc,npc),temp,
     & sn(*)
      real var(nt),tf_se(nt,nt),varsum
      integer i1(n1),i2(n2)

ccc   convert b into array TF (using PCs, var)
      do i=1,n2
         do j = 1,n2
            tf(i2(i),i2(j)) = (0.,0.)
         enddo
         do j = 1,n1
            temp = (0.,0.)
            do k = 1,npc
               temp = temp + conjg(u(j,k))*b(k,i)
            enddo
            tf(i2(i),i1(j)) = temp/sqrt(var(i1(j))+varsum)
         enddo
      enddo

ccc   now TF variances
      do i=1,n2
         do j = 1,n2
            tf_se(i2(i),i2(j)) = 0.
         enddo
         ii = (i*(i+1))/2
         do j = 1,n1
            temp = (0.0,0.0)
            do k = 1,npc
               do l = 1,npc
                  temp = temp + conjg(u(j,k))*xxinv(k,l)*u(j,l)
               enddo
            enddo
            tf_se(i2(i),i1(j)) = sqrt(real(temp)*sn(ii)
     &                                  /(var(i1(j))+varsum))
         enddo
      enddo

      return
      end
ccc_____________________________________________________________________
ccc
      function ipt_sn(nch_grp,ista)
      integer ipt_sn,nch_grp(*),ista,itemp
      itemp  = 1
      do k = 1,ista-1
         itemp = itemp + (nch_grp(k)*(nch_grp(k)+1))/2
      enddo
      ipt_sn = itemp 
      return
      end
