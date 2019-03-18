      subroutine cdrlud(a,b,n,mx,m,convrg)
c----------------------------------------------------------------------
c     Given a complex*16 matrix a, this returns its inverse matrix.
c     After calling of this program a is replaced by the lu
c     decomposition. The half pivoting is employed.
c     coded by H.Akai, 1 Dec. 1996, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (nmax=500)
      complex*16 a(n,n,mx),b(n,n,mx)
      integer indx(nmax)
      logical convrg(mx)
      if(n .gt. nmax) call errtrp(1,'cdrlud','n too large')
      call sbtime(2,0)
      do 40 k=1,m
      if(.not. convrg(k)) then
      do 10 i=1,n
      do 20 j=1,n
   20 b(i,j,k)=(0d0,0d0)
   10 b(i,i,k)=(1d0,0d0)
      call cludcm(a(1,1,k),n,indx,d)
      do 30 j=1,n
   30 call clubks(a(1,1,k),n,indx,b(1,j,k))
      a(1,1,k)=a(1,1,k)*d
      endif
   40 continue
      call sbtime(2,1)
      end
