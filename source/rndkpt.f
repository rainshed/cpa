      subroutine rndkpt(korder,nk,iwk)
c---------------------------------------------------------------------
c     Rearragne k-points in the BZ such that they distribute
c     uniformly in Weyl's criterion in every step of sampling.
c     Coded H.Akai, 13 Sep. 1996, Osaka
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer korder(nk),iwk(nk)
      data irnd/0/
      if(irnd .eq. 0) then
      do 30 k=1,nk
   30 korder(k)=k
      return
      else
      call errtrp(3,'rndkpt','randomization takes place')
      c=sqrt(3d0)
      do 10 k=1,nk
   10 iwk(k)=0
      i=0
      do 20 n=1,1000000
      if(i .ge. nk) return
      x=c*dble(n)
      ix=x
      k=(x-dble(ix))*dble(nk)+1d0
      k=min(k,nk)
      if(iwk(k) .eq. 0) then
      i=i+1
      iwk(k)=i
      korder(i)=k
      endif
   20 continue
      endif
c      call errtrp(1,'rndkpt','nk too large')
      end
