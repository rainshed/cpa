      function cpolin(x,xin,yin,n,err)
c----------------------------------------------------------------------
c     Given arrays xin and yin, each of length n and given value x,
c     this function returns a value polin0. If p(x) is the polynomial
c     of degree n-1 such that p(xin(i))=yin(i), i=1,..,n, then the
c     program returned value polin0=p(x). Complex*16 version.
c     See 'Numerical Recipes', chapter 3, for further details.
c     coded by H.Akai, Jan. 1991, Osaka
c     modified by H. Akai, Tokyo, Jan. 2018
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 xin(n),yin(n),x,cd,cpolin
      complex*16,allocatable::c(:),d(:)
      data small/1d-30/
      allocate(c(n),d(n))
      do 10 i=1,n
      c(i)=yin(i)
   10 d(i)=yin(i)
      cpolin=yin(1)
      do 20 j=1,n-1
      do 30 i=1,n-j
c     ---an error occurs if two xin's come too closer.
      if(abs(xin(i)-xin(i+j)) .lt. small)
     &   call errtrp(1,'cpolin','data error')
c     ---recursive relation eq.(3.1.5) of Numerical Recipes.
      cd=(c(i+1)-d(i))/(xin(i)-xin(i+j))
      d(i)=(xin(i+j)-x)*cd
   30 c(i)=(xin(i)-x)*cd
   20 cpolin=cpolin+c(1)
      err=min(abs(c(1)),abs(d(1)))
      deallocate(c,d)
      end
