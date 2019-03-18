      subroutine cgnwt(e,wt,ng,kmx,ew,ez)
c-----------------------------------------------------------------------
c     Gives weighting function which is used for the contour
c     integration along complex energy.
c     coded by H.Akai, 1983, Juelich
c     Modified by H. Akai, Tokyo, Jan. 2018.
c-----------------------------------------------------------------------
      parameter(ngmx=21)
      implicit real*8 (a-h,o-z)
      complex*16 e(kmx),wt(ng,3,kmx),t(ngmx),t0(ngmx),t1(ngmx),t2(ngmx)
      if(ng .gt. ngmx) call errtrp(1,'cgnwt','ng too large')
      do 20 k=1,kmx
      t(1)=(1d0,0d0)
      t(2)=(e(k)-ew)/ez
      t0(1)=t(2)
      t0(2)=t(2)**2/2d0
      t1(1)=t0(2)
      t1(2)=t(2)**3/3d0
      do 30 i=3,ng
      t(i)=2d0*t(2)*t(i-1)-t(i-2)
      t0(i)=(dble(i-4)*t0(i-2)-(2d0-4d0*t0(2))*t(i-1))/dble(i)
   30 t1(i)=(dble(i-1)*t0(i-1)-(1d0-2d0*t0(2))*t(i))/dble(i+1)
      t2(1)=t1(2)
      do 40 i=2,ng
   40 t2(i)=(dble(i-1)*t1(i-1)+t0(i)-t(2)*(1d0-2d0*t0(2))*t(i))
     &       /dble(i+2)
      do 20 i=1,ng
      wt(i,1,k)=ez*t0(i)
      wt(i,2,k)=ez**2*t1(i)+ew*wt(i,1,k)
   20 wt(i,3,k)=ez**3*t2(i)+2d0*ew*ez**2*t1(i)+ew**2*wt(i,1,k)
      end
