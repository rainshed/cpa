      subroutine poisna(r,v,anclr,dr,xr,meshr)
c-----------------------------------------------------------------------
c     This subroutine solves Poisson equation by use of the
c     Simpson's integration method. It is rather straightfoward
c     and not an elegant way. The advantage of this method is that
c     it is quite tough and hardly failes even in the case of
c     singular charge distribution such as muonic atoms.
c     coded by H.Akai, 1988, Osaka
c     revised by H.Akai, 1992, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 r(meshr),v(meshr),dr(meshr),xr(meshr),old,mid,new
c     ---i first integrate the contribution of outside.
      v(meshr-1)=0d0
      old=xr(meshr-1)*dr(meshr-1)*r(meshr-1)
      do 10 k=meshr-3,1,-2
      mid=xr(k+1)*dr(k+1)*r(k+1)
      new=xr(k)*dr(k)*r(k)
      v(k+1)=v(k+2)+(5d0*old-new+8d0*mid)/6d0
      v(k)  =v(k+1)+(5d0*new-old+8d0*mid)/6d0
   10 old=new
c     ---then i integrate the contribution of inside.
      z=2d0*(xr(1)**3*r(1)/3d0-anclr)
      v(1)=v(1)+z/xr(1)
      old=xr(1)**2*dr(1)*r(1)
      do 20 k=1,meshr-3,2
      mid=xr(k+1)**2*dr(k+1)*r(k+1)
      new=xr(k+2)**2*dr(k+2)*r(k+2)
      z=z+(5d0*old-new+8d0*mid)/6d0
      v(k+1)=v(k+1)+z/xr(k+1)
      z=z+(5d0*new-old+8d0*mid)/6d0
      v(k+2)=v(k+2)+z/xr(k+2)
   20 old=new
      return
      end
