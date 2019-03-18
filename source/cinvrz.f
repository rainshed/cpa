      subroutine cinvrz(p,q,n,mx,m)
c----------------------------------------------------------------------
c     Given a complex matrix p, this returns its inverse matrix.
c     After calling of this program p is replaced by the lu
c     decomposition. Half pivoting is employed.
c     In this version the operations running over the first index of
c     the matrix come up to the outermost loop, suitable for a big
c     matrix with rather small size regardint the first index.
c     coded by H.Akai., Feb. 4, 1992, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(np=90)
      complex*16 p(mx,n,n),q(mx,n,n),a(np,np),b(np,np)
      integer indx(np)
      if(n .gt.np) call errtrp(1,'cinvrz','n too large')
      do 40 k=1,m
      do 10 i=1,n
      do 20 j=1,n
      a(i,j)=p(k,i,j)
   20 b(i,j)=(0d0,0d0)
   10 b(i,i)=(1d0,0d0)
      call ludcmp(a,n,np,indx,d)
      do 30 j=1,n
   30 call lubksb(a,n,np,indx,b(1,j))
      a(1,1)=a(1,1)*d
      do 40 i=1,n
      do 40 j=1,n
      p(k,i,j)=a(i,j)
   40 q(k,i,j)=b(i,j)
      end
