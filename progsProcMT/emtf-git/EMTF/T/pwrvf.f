       subroutine pwrvf(x,nch,ns,nf,p,wrk)
       complex x(nch,ns,nf)
       real p(nf),wrk(nch,0:nf)

c************given the array x (nch channels , ns time segments, nf
c                   frequencies )  find average power for each
c                frequency, and return a damped inverse of this
c                crude average power spectrum for whitening the
c                 spectrum

c........find average power for each freq, channel
       do 10 i = 1,nch
       do 10 j = 1,nf
       wrk(i,j) = 0.
          do 5 k = 1,ns
5         wrk(i,j) =wrk(i,j) + x(i,k,j)*conjg(x(i,k,j))
10     wrk(i,j) = wrk(i,j)/ns

c.........normalize so that each channel has avg pwr 1 (avgd over freq)
       do 20 i = 1,nch
       wrk(i,0) = 0.
         do 15 j = 1,nf
15       wrk(i,0) = wrk(i,0)+wrk(i,j)
       wrk(i,0) = wrk(i,0)/nf
         do 18 j = 1,nf
         wrk(i,j) = wrk(i,j)/wrk(i,0)
18       continue
20    continue

c......... find average power for each freq.
      do 30 j = 1,nf
      p(j) =0.
         do 25 i = 1,nch
25       p(j) = p(j) + wrk(i,j)
      p(j) = p(j)/nch
      wrk(1,j) = p(j)
30    continue

c.......... find median average power
      call sort(p,nf)

c.......... add to average power to determine relative weights
      pmed = p(nf/2)
      do 40 j = 1,nf
40    p(j) = 1./(wrk(1,j) + pmed)
      return
      end
