      subroutine cubmtr(u,mxl)
c-----------------------------------------------------------------------
c     Generate cubic rotation matrices for real harmonics
c     Coded by H. Akai, Dec 96, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 u((2*mxl-1)**2*mxl,24)
      pi=4d0*atan(1d0)
      hpi=5d-1*pi
      hpin=-hpi
c
c---  E identity operation
      call getrtm(u(1,1),0d0,0d0,0d0,mxl)
c
c---  C_4 (1/2 pi rotation around x, y, and z axes)
      call getrtm(u(1,2),hpi,hpi,hpin,mxl)
      call getrtm(u(1,3),0d0,hpi,0d0,mxl)
      call getrtm(u(1,4),hpi,0d0,0d0,mxl)
c
c---  C_2 (pi rotation around x, y, and z axes)
      call prdrtm(u(1,2),u(1,2),u(1,5),mxl)
      call prdrtm(u(1,3),u(1,3),u(1,6),mxl)
      call prdrtm(u(1,4),u(1,4),u(1,7),mxl)
c
c---  C_4^3 (3/2 pi rotation around x, y, and z axes)
      call prdrtm(u(1,2),u(1,5),u(1,8),mxl)
      call prdrtm(u(1,3),u(1,6),u(1,9),mxl)
      call prdrtm(u(1,4),u(1,7),u(1,10),mxl)
c
c---  C_3 (2/3 pi rotation around 111, 11-1, 1-11, and 1-1-1 axes)
      call getrtm(u(1,11),hpi,hpi,0d0,mxl)
      call getrtm(u(1,12),0d0,hpi,hpin,mxl)
      call getrtm(u(1,13),0d0,hpin,hpi,mxl)
      call getrtm(u(1,14),hpin,hpin,0d0,mxl)
c
c---  C_3^2 (4/3 pi rotation around 111, 11-1, 1-11, and 1-1-1 axes)
      call prdrtm(u(1,11),u(1,11),u(1,15),mxl)
      call prdrtm(u(1,12),u(1,12),u(1,16),mxl)
      call prdrtm(u(1,13),u(1,13),u(1,17),mxl)
      call prdrtm(u(1,14),u(1,14),u(1,18),mxl)
c
c---  C_2' (pi rotation around 110, 1-10, 011, 01-1, 101, 10-1 axes)
      call getrtm(u(1,19),hpi,pi,0d0,mxl)
      call getrtm(u(1,20),hpin,pi,0d0,mxl)
      call getrtm(u(1,21),hpi,hpi,hpi,mxl)
      call getrtm(u(1,22),hpi,hpin,hpi,mxl)
      call getrtm(u(1,23),0d0,hpin,pi,mxl)
      call getrtm(u(1,24),0d0,hpi,pi,mxl)
c
      end
