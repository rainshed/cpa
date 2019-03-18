      function polin0(x,xa,ya,n,dy)
c----------------------------------------------------------------------
c     Given arrays xa and ya, each of length n and given value x,
c     this function returns a value polin0. if p(x) is the polynominal
c     of degree n-1 such that p(xa(i))=ya(i), i=1,..,n, then
c     the returned value polin0=p(x).
c     See 'Numerical Recipes', chapter 3
c     Coded by H. Akai, Jan. 1991, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (nmx=20)
      dimension xa(n),ya(n),c(nmx),d(nmx)
      data small/1d-30/
      ns=1
      dif=abs(x-xa(1))
      do 40 i=1,n
      dift=abs(x-xa(i))
      if(dift .lt. dif) then
      ns=i
      dif=dift
      endif
      c(i)=ya(i)
   40 d(i)=ya(i)
      polin0=ya(ns)
      ns=ns-1
      do 60 m=1,n-1
      do 50 i=1,n-m
      ho=xa(i)-x
      hp=xa(i+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
c
c     an error can occur if two xa's are identical
      if(abs(den) .lt. small) then
      call errtrp(2,'polin0','data error')
      return
      endif
c
      den=w/den
      d(i)=hp*den
   50 c(i)=ho*den
      if(2*ns .lt. n-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
   60 polin0=polin0+dy
      return
      end
