c
c************************************
c
      subroutine pterst(irec1,ifirst,ifirstd,next,lstart,lstartd)

c    initializes pointers for decimation routines

      include 'decimate.inc'
      integer ifirst(nd),ifirstd(nd),next(nd),nmod(10)
      logical lstart(nd),lstartd(nd)

        nmod(1)=1
        do 100 id = 1,nd

c       set up; set pointer next

        lstart(id) = .false.
        lstartd(id) = .false.
        next(id) = 0
        r = nmod(id)*(nwin(id)-olap(id))

c       find first starting record for level id; set pointer ifirst

        i1 = irec1
        i2 = i1 + nint(r) + 1
        i3 = nmod(id)
        ii = 0
           do 20 i = i1,i2,i3
c  changed according to message from Gary Egbert 12/16/97           
c          k = nint((mod(float(i),r)-ioffs(id))/nmod(id))
           k = nint((mod(dfloat(i),dble(r))-ioffs(id))/nmod(id))
           if (k .eq. 0) go to 25
20         ii = ii + 1

25      ifirst(id) = ii

        if(id.lt.nd) then

c       find first record for level id+1; set decimation pointer ifirstd

           nmod(id+1) = nmod(id)*idec(id+1)
           i1 = irec1 + (nfc(id+1)-1)*nmod(id)
           i2 = i1 + nmod(id+1)
           i3 = nmod(id)
           ii = nfc(id+1) - 1

              do 30 i = i1,i2,i3
              if(mod(i,nmod(id+1)).eq. ioffr(id+1)) go to 35
30            ii = ii + 1
35         irec1 = i
           ifirstd(id) = ii
           end if

100     continue
        return
        end
