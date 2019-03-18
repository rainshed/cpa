      subroutine prdrtm(u1,u2,u3,mxl)
c-----------------------------------------------------------------------
c     make the product of the rotaion matrices u1 and u2 and put the
c     result into u3.
c     coded by H.Akai, 14 Dec. 1996, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 u1((2*mxl-1)**2,mxl),u2((2*mxl-1)**2,mxl)
     &      ,u3((2*mxl-1)**2,mxl)
      do 10 l=1,mxl
      n=2*l-1
   10 call prdmtr(u1(1,l),u2(1,l),u3(1,l),n)
      end
