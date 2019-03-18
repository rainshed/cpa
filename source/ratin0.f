      function ratin0(x,xa,ya,n,dy)
c----------------------------------------------------------------------
c     Given arrays xa and ya, each of length n and given value x,
c     this function returns a value ratin0. the value returned is that
c     of the diagonal rational function, evaluated at x, which passes
c     through the n points (xa(i),y(i)),i=1,..,n.
c     See Numerical Recipes pp83-85.
c     coded by H.Akai, Feb. 7, 1991, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (nmx=30)
      dimension xa(n),ya(n),c(nmx),d(nmx)
      data small/1d-30/
      if(n .gt. nmx) call errtrp(1,'ratin0','n too large')
      ns=1
      hh=abs(x-xa(1))
      do 40 i=1,n
      h=abs(x-xa(i))
      if(h .lt. small) then
      ratin0=ya(i)
      dy=0d0
      return
      else if(h .lt. hh) then
      ns=i
      hh=h
      endif
      c(i)=ya(i)
   40 d(i)=ya(i)+small
      ratin0=ya(ns)
      ns=ns-1
      do 60 m=1,n-1
      do 50 i=1,n-m
      w=c(i+1)-d(i)
      h=xa(i+m)-x
      t=(xa(i)-x)*d(i)/h
      dd=t-c(i+1)
c
c     an error can occur if the interpolating function has a pole
c     at the requested value of x.
      if(abs(dd) .lt. small) call errtrp(1,'ratin0','pole at x')
c
      dd=w/dd
      d(i)=c(i+1)*dd
   50 c(i)=t*dd
      if(2*ns .lt. n-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
   60 ratin0=ratin0+dy
      return
      end
