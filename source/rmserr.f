      subroutine rmserr(v1,v2,rms,dr,xr,meshr)
c-----------------------------------------------------------------------
c     rms error analysis
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension dr(meshr),xr(meshr),v1(meshr),v2(meshr)
      atvol=xr(meshr)**3/3d0
      volins=xr(meshr-1)**3/3d0
      volint=atvol-volins
      rms=0d0
      ints=meshr-2
c     ker=0
c     chk=0d0
c     do 30 k=1,meshr-1
c     dif=((v1(k)-v2(k))*xr(k))**2
c     if(dif .gt. chk) then
c     chk=dif
c     ker=k
c     endif
c  30 continue
c     if(ker .ne. 0)
c    &  write(*,'(1x,a,i3,1p,e13.6)')'   max err at',ker,chk
      do 10 k=2,ints,2
   10 rms=rms+((v1(k)-v2(k))*xr(k))**2*dr(k)
      rms=rms*2d0
      do 20 k=3,ints,2
   20 rms=rms+((v1(k)-v2(k))*xr(k))**2*dr(k)
      rms=rms*2d0+((v1(1)-v2(1))*xr(1))**2*dr(1)
     &         +((v1(meshr-1)-v2(meshr-1))*xr(meshr-1))**2
     &         *dr(meshr-1)
      rms=rms/3d0+volint
     &         *(v1(meshr)-v2(meshr))**2
c     --- Here, rms is the root-mean-square error weighted by
c         the volume. Therefore, this quantity will proportional
c         to the volume over which the root-mean-square is
c         associated.
      rms=atvol*sqrt(rms/atvol)
      return
      end
