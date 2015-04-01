      subroutine n_rbst(xr,xc,xp,nd,nt,n_grp,ih_grp,nch_grp,
     &   sndiag,s,period,itmax_stk,itmax_ln,itmax_cln,r0,ioux,
     &   iocm,iotf)
   
ccc    NEW ROBUST VERSION OF INCOHERENT NOISE VARIANCE ESTIMATION
ccc      ROUTINE
ccc     estimates of "noise" in each channel, using methods of 
ccc     Egbert, 1996  ;  Channels are divided into n_grp groups;
ccc      group i contatins nch_grp(i) consecutive channels beginnig
ccc      with channel ih_grp(i)
ccc     Noise estimates are based on fitting all channels in a group
ccc      to PCs constructed from all channels in the other groups.
ccc     After computation of raw residual variances, a linear system
ccc       is solved (in the call to var_adj) to compute approximately
ccc       unbiased estimates of incoherent noise variances for each
ccc       channel

ccc      VARIOUS FLAGS WHICH CONTROL OPTIONS:
ccc      iocm is output unit for correlation matrix output (if > 0 ...
ccc        otherwise this is not output)
ccc      iotf gives output unit for prediction TF (if > 0 ...
ccc        otherwise this is not output)
ccc      ioux gives output unit for residuals (if >0 )

      include 'iosize.inc'
      parameter (ns1mx=(nchmx*(nchmx+1))/2,nwork=ntmx*ntmx*3)
      parameter (ntmx_sq = ntmx*ntmx)

      integer ih_grp(ntmx),nch_grp(ntmx),n_grp
      complex xr(nd,nt),xc(nd,nt),xp(nd,nt),s(*)
      real sndiag(ntmx),tf_se(ntmx,ntmx)
      complex sn(ntmx_sq),tf(ntmx,ntmx),work(nwork)
      integer iotf
      character*80 comment

ccc   imxc1 is # of iterations for robust regression on PCs estimate
      imxc1 = itmax_cln
      do i = 1,nt
         sndiag(i) = 0.0
      enddo

ccc   copy raw data (xr) into cleaned data array (xc)
      call datawt(0,nd,nt,xr,xc,xc,sn,r0)
      if(ioux.gt.0) then
ccc      output raw input data
         write(comment,*) 'Raw data in n_rbst'
         call wrtx(ioux,xc,nt,nd,1,sndiag,comment)
      endif

ccc   initial estimate of noise in each data channel
      call ln_rbst(xr,xc,xp,sndiag,s,nt,nd,n_grp,ih_grp,sn,tf,tf_se,
     &  nch_grp,itmax_stk,itmax_ln)

ccc   copy raw data (xr) into cleaned data array (xc)
      call datawt(0,nd,nt,xr,xc,xc,sn,r0)
ccc   extract diagonal of residual covariance estiamte
c       write(*,*) 'n_grp',n_grp
c       write(*,*) 'ih_grp',(ih_grp(i),i=1,n)
      do i=1,n_grp
         isn = ipt_sn(nch_grp,i)
         call diag(sn(isn),nch_grp(i),sndiag(ih_grp(i)))
      enddo
c      write(*,*) 'sndiag',(sndiag(l),l=1,14)

ccc   iterative refinement   
      do iter=1,imxc1
         call ln_rbst(xr,xc,xp,sndiag,s,nt,nd,n_grp,ih_grp,sn,tf,tf_se,
     &     nch_grp,itmax_stk,itmax_ln)
         do i=1,n_grp
            isn = ipt_sn(nch_grp,i)
            call diag(sn(isn),nch_grp(i),sndiag(ih_grp(i)))
         enddo
c         write(*,*) 'iter = ',iter
c         write(*,*) 'sndiag',(sndiag(l),l=1,14)
ccc      using raw data (XR) + predictions (XP) form cleaned data (XC)
         call datawt(iter,nd,nt,xr,xp,xc,sndiag,r0)
      enddo ! iter

ccc   optional output of estimated TF
      if(iotf.gt.0) call wrt_tf(iotf,period,tf,tf_se,nt,nd)

ccc   variance adjustment (see Appendix A in Egbert, 1996)
      call var_adj(nt,tf,sndiag,work)
         
ccc   recompute SDM (no iterations, use cleanded data...
ccc           ... appropriate for comparison with noise estimates
      call rbstk2(xc,nt,nd,s,xp,0)

ccc   optional output of correlation matrix, group residual correlation
ccc   matrices
      if(iocm .gt. 0) then
        write(iocm,*) '************  PERIOD = ',period,'  ************'
        write(iocm,*) '************* # of Data = ',nd,'   *************'
        call wrt_cormat(iocm,s,nt)
        call wrt_sn(iocm,sn,nch_grp,n_grp)
      endif
      return
      end
c______________________________________________________________________
c
      subroutine wrt_tf(iotf,period,tf,tf_se,nt,nd)
      complex tf(nt,nt)
      real tf_se(nt,nt),period
      integer iotf
      write(iotf,'(e12.4,i12)') 1./period,nd
      do i = 1,nt
         do j=1,nt
            write(iotf,'(3e12.4)') tf(i,j),tf_se(i,j)
         enddo
      enddo
      return
      end
c_____________________________________________________________________
c
      subroutine mkgrp(option_grp,ih,ie,nsta,nt,ih_grp,n_grp,nch_grp)

ccc   makes arrays ih_grp,n_grp,nch_grp for call to n_rbst based
ccc   on value of character variable 

      character*1 option_grp
      integer ih(*),ie(*),nsta,nt,ih_grp(*),n_grp,nch_grp(*)

      if(option_grp.eq.'S') then
ccc      groups are all channels at each site
         n_grp = nsta
         do ista = 1,nsta
            nch_grp(ista) = ih(ista+1) - ih(ista)
            ih_grp(ista) = ih(ista)
         enddo
         ih_grp(nsta+1) = nt+1
      else if(option_grp.eq.'A') then
ccc      each channel is a group by itself
         n_grp = nt
         do i = 1,nt+1
            nch_grp(i) = 1
            ih_grp(i) = i
          enddo
      else
ccc      use all channels of a component type at a fixed site as a group
ccc      (default option)
         ii = 1
         do ista = 1,nsta
            nch_grp(ii) = ie(ista)-ih(ista)
            if(nch_grp(ii).gt.0) then
               ih_grp(ii) = ih(ista)
               ii = ii + 1
             endif
             nch_grp(ii) = ih(ista+1)-ie(ista)
             if(nch_grp(ii).gt.0) then
                ih_grp(ii) = ie(ista)
                ii = ii + 1
             endif
          enddo
          ih_grp(ii) = nt+1
          n_grp= ii - 1
       endif
       return
       end
