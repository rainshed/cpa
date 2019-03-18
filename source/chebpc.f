      subroutine chebpc(c,d,n)
c-----------------------------------------------------------------------
c  (C) Copr. 1986-92 Numerical Recipes Software !+!).
c     modified and addapted to KKR packege by H.Akai, 1993
c     modified by H. Akai, Tokyo, Jan. 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 c(n),d(n)
      real*8,allocatable::dd(:)
      allocate(dd(n))
      do 10 j=1,n
      d(j)=0d0
   10 dd(j)=0d0
      d(1)=c(n)
      do 20 j=n-1,2,-1
      do 30 k=n-j+1,2,-1
      sv=d(k)
      d(k)=2d0*d(k-1)-dd(k)
   30 dd(k)=sv
      sv=d(1)
      d(1)=-dd(1)+c(j)
   20 dd(1)=sv
      do 40 j=n,2,-1
   40 d(j)=d(j-1)-dd(j)
      d(1)=-dd(1)+c(1)
      deallocate(dd)
      end
