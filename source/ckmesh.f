      subroutine ckmesh(rtin,ratm,xr,dr,vu,vd,wk,meshr)
c----------------------------------------------------------------------
c     This program check radial mesh read in from file if it is
c     consistent with the lattice parameter given previously. If
c     not, it will interpolate the potential data such that they
c     are compatible with the present lattice parameters.
c     coded by H.Akai, 1986, Juelich
c     revised by H.Akai, April 1992, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 xr(meshr),dr(meshr),wk(meshr,3),vu(meshr),vd(meshr)
      data small/1d-6/
      call rmesha(1d-6,rtin,ratm,wk(1,2),wk(1,1),meshr)
      if(abs(ratm-xr(meshr)) .gt. small .or.
     &         abs(rtin-xr(meshr-1)) .gt. small) then
c     write(6,1000)
c1000 format('   ***msg ckmesh...new mesh generated')
      do 10 k=1,meshr-1
   10 wk(k,3)=vu(k)*xr(k)
      do 20 k=1,meshr-1
   20 vu(k)=polint(wk(k,1),xr,wk(1,3),meshr-1)/wk(k,1)
      do 30 k=1,meshr-1
   30 wk(k,3)=vd(k)*xr(k)
      do 40 k=1,meshr-1
   40 vd(k)=polint(wk(k,1),xr,wk(1,3),meshr-1)/wk(k,1)
      endif
c     ---anyway, the new mesh consistent with rmt and ratm is generated
c        to avoid possible numerical errors.
      call equarr(wk(1,1),xr,meshr)
      call equarr(wk(1,2),dr,meshr)
      return
      end
