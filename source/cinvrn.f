      subroutine cinvrn(a,b,nx,n,m)
c----------------------------------------------------------------------
c     Given a complex matrix a, this returns its inverse matrix.
c     After calling of this program a is replaced by the lu
c     decomposition. No pivoting is employed.
c     In this version the operations running over the first index of
c     the matrix come up to the outermost loop, suitable for a big
c     matrix with rather small size regardint the first index.
c     coded by H.Akai, Feb. 4, 1992, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(nx,nx,m),b(nx,nx,m),c
      do 10 k=1,m
      do 20 i=1,n-1
      c=-1d0/a(i,i,k)
      do 20 j=i+1,n
      a(j,i,k)=c*a(j,i,k)
      do 30 l=1,i-1
   30 a(j,l,k)=a(j,l,k)+a(j,i,k)*a(i,l,k)
      do 20 l=i+1,n
   20 a(j,l,k)=a(j,l,k)+a(j,i,k)*a(i,l,k)
      do 10 i=n,1,-1
      c=1d0/a(i,i,k)
      do 40 j=1,i-1
      b(i,j,k)=a(i,j,k)
      do 40 l=i+1,n
   40 b(i,j,k)=b(i,j,k)-a(i,l,k)*b(l,j,k)
      do 50 j=i+1,n
      b(i,j,k)=(0d0,0d0)
      do 50 l=i+1,n
   50 b(i,j,k)=b(i,j,k)-a(i,l,k)*b(l,j,k)
      b(i,i,k)=(1d0,0d0)
      do 60 l=i+1,n
   60 b(i,i,k)=b(i,i,k)-a(i,l,k)*b(l,i,k)
      do 10 j=1,n
   10 b(i,j,k)=c*b(i,j,k)
      end
