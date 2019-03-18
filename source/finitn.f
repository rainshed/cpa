      subroutine finitn(v,anclr,xr,meshr)
c-----------------------------------------------------------------------
c     For given anclr, this program makes a finite nuclear size
c     correction on the potential.
c     Coded by H. Akai, Nov 1985, Osaka
c     Revised by H. Akai, 27 May, 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(meshr),xr(meshr)
c     return
      a=atmmas(anclr)
      if(abs(a) .lt. 1d-20) return
      r=.226767d-4*a**(1d0/3d0)
c     write(*,'(1x,a,1p,3e13.6)')'r=',r
      do 10 k=1,meshr-1
      if(xr(k) .gt. r) return
   10 v(k)=v(k)+2d0*anclr*(r/xr(k)+5d-1*(xr(k)/r)**2-1.5d0)/r
c  10 write(*,'(1x,i3,1p,3e13.6)')k,xr(k),v0*xr(k),v(k)
      end
