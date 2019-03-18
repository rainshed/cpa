      subroutine tchmta(ew,ez,elvl,tm,ng)
c-----------------------------------------------------------------------
c     Generate Thcebycheff transformation matrix.
c     coded by H.Akai, 1979, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension elvl(ng),tm(ng,ng)
      data pai/3.14159265358d0/
      p1=pai/(2d0*dble(ng))
      do 20 k=1,ng
      ang=dble(2*k-1)*p1
      tm(1,k)=1d0
      tm(2,k)=cos(ang)
      elvl(k)=ew+ez*tm(2,k)
      do 10 l=3,ng
   10 tm(l,k)=2d0*tm(2,k)*tm(l-1,k)-tm(l-2,k)
      tm(1,k)=tm(1,k)/dble(ng)
      do 20 l=2,ng
   20 tm(l,k)=2d0*tm(l,k)/dble(ng)
      return
      end
