      subroutine cinvrx(a,b,n,mx,m,ibsadr,kmx,convrg)
c----------------------------------------------------------------------
c     Given a complex matrix a, this returns its inverse matrix.
c     After calling of this program a is replaced by the lu
c     decomposition. No pivoting is employed.
c     In this version the operations running over the first index of
c     the matrix come up to the outermost loop, suitable for a big
c     matrix with rather small size regardint the first index.
c     coded by H.Akai., Feb. 4, 1992, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(n,n,mx),b(n,n,mx)
      logical convrg(kmx)
      call sbtime(2,0)
      do 40 k=1,m
      iadr=ibsadr+k
      kk=mod(iadr-1,kmx)+1
      if(.not. convrg(kk)) then
      do 10 i=1,n-1
      do 10 j=i+1,n
      a(j,i,k)=-a(j,i,k)/a(i,i,k)
      do 30 l=1,i-1
   30 a(j,l,k)=a(j,l,k)+a(j,i,k)*a(i,l,k)
      do 10 l=i+1,n
   10 a(j,l,k)=a(j,l,k)+a(j,i,k)*a(i,l,k)
      do 80 i=n,1,-1
      do 60 j=1,i-1
      b(i,j,k)=a(i,j,k)
      do 60 l=i+1,n
   60 b(i,j,k)=b(i,j,k)-a(i,l,k)*b(l,j,k)
      do 70 j=i+1,n
      b(i,j,k)=(0d0,0d0)
      do 70 l=i+1,n
   70 b(i,j,k)=b(i,j,k)-a(i,l,k)*b(l,j,k)
      b(i,i,k)=(1d0,0d0)
      do 90 l=i+1,n
   90 b(i,i,k)=b(i,i,k)-a(i,l,k)*b(l,i,k)
      do 80 j=1,n
   80 b(i,j,k)=b(i,j,k)/a(i,i,k)
      endif
   40 continue
      call sbtime(2,1)
      return
      end
