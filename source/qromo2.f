      subroutine qromo2(da,e,l,a,b,ss)
c-----------------------------------------------------------------------
c     Romberg integration on an open interval.
c     See Numerical Recipes 2nd edition p135.
c     Modified by H. Akai, 13 April 1996, osaka
c  (C) Copr. 1986-92 Numerical Recipes Software !+!).
c      Modified to adapt OpenMP version by H. Akai, 6 Jan 2016, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (eps=1d-6, jmax=14, jmaxp=jmax+1, k=5, small=1d-8)
      real*8 h(jmaxp),s(jmaxp)
      km=k-1
      h(1)=1d0
      do 10 j=1,jmax
      call mdexp2(da,e,l,a,b,s(j),j)
      if(j .ge.3 .and. abs(s(j)) .lt. small) then
      ss=0d0
      return
      endif
      if (j .ge. k) then
      ss=polin0(0d0,h(j-km),s(j-km),k,dss)
      if (abs(dss) .le. eps*abs(ss)) return
      endif
      s(j+1)=s(j)
      h(j+1)=h(j)/9d0
   10 continue
      call errtrp(2,'qromo2','too many steps')
      end
