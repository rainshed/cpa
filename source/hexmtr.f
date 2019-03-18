      subroutine hexmtr(u,mxl)
c-----------------------------------------------------------------------
c     Generate hexagonal rotation matrises for real harmonics
c     Coded by H. Akai, Dec 96, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 u((2*mxl-1)**2*mxl,24)
      pi=4d0*atan(1d0)
      hpi=5d-1*pi
      hpin=-hpi
      tpi=pi/3d0
      tpin=-pi/3d0
      spi=pi/6d0
      spin=-spi
c
c---  E identity operation
      call getrtm(u(1,1),0d0,0d0,0d0,mxl)
c
c---  C_6 (1/3 pi rotation around z axis)
      call getrtm(u(1,2),tpi,0d0,0d0,mxl)
c---  C_3=C_6^2
      call prdrtm(u(1,2),u(1,2),u(1,3),mxl)
c---  C_2=C_6^3
      call getrtm(u(1,4),pi,0d0,0d0,mxl)
c---  C_3-=C_6^4
      call prdrtm(u(1,4),u(1,2),u(1,5),mxl)
c---  C_6-=c_6^5
      call prdrtm(u(1,4),u(1,3),u(1,6),mxl)
c---  C_2' (pi rotaion around y=x tan(n pi/6) (n=0,5))
      call getrtm(u(1,7),hpi,pi,hpin,mxl)
      call getrtm(u(1,8),tpi,pi,tpin,mxl)
      call getrtm(u(1,9),spi,pi,spin,mxl)
      call getrtm(u(1,10),0d0,pi,0d0,mxl)
      call getrtm(u(1,11),spin,pi,spi,mxl)
      call getrtm(u(1,12),tpin,pi,tpi,mxl)
      call clrarr(u(1,13),12*(2*mxl-1)**2*mxl)
c
      end
