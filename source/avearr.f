      subroutine avearr(a,b,n)
c-----------------------------------------------------------------------
c     give spin up/down averaged data
c     coded by H.Akai, 1993, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      do 10 k=1,n
      a(k)=5d-1*(a(k)+b(k))
   10 b(k)=a(k)
      return
      end
