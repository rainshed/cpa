      subroutine cprinv(a,nx,n,m)
c----------------------------------------------------------------------
c     Given a complex matrix a, this returns its lu-decomposition.
c     No pivoting is used.
c     Coded by H.Akai., Feb. 4, 1992, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(nx,nx,m)
      do 30 k=1,m
      do 10 i=1,n-1
      do 10 j=i+1,n
      a(j,i,k)=-a(j,i,k)/a(i,i,k)
      do 20 l=1,i-1
   20 a(j,l,k)=a(j,l,k)+a(j,i,k)*a(i,l,k)
      do 10 l=i+1,n
   10 a(j,l,k)=a(j,l,k)+a(j,i,k)*a(i,l,k)
   30 continue
      end
