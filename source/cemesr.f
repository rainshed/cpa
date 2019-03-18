      subroutine cemesr(ef,ewidth,edelt,ebtm,e,kmx)
c----------------------------------------------------------------------
c     Generate a energy contour along the real axis.
c     The following are the examples of some typical cases.
c       ewidth/2.5d0/, edelt/3d-3/, ref/1.0d0/
c       ewidth/0.500d0/, edelt/1d-4/, ref/0.9d0/
c       ewidth/1.500d0/, edelt/1d-2/, ref/0d0/
c       ewidth/3.300d0/, edelt/3d-3/, ref/0.9d0/
c     coded by H.Akai, 1983, Juelich
c     latest version 30 Nov. 1997, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension e(kmx)
      complex*16 e,di
c     data ref/0.9d0/
      data ref/0.75d0/
c     data ref/1d0/
c     data ref/0.666666d0/
c     data ref/0.5d0/
c     data ref/0d0/
      di=dcmplx(0d0,edelt)
      de=ewidth/dble(kmx-1)
      kef=dble(kmx-1)*ref+1
      el=ef-de*dble(kef-1)
      do 10 k=1,kmx
   10 e(k)=el+dble(k-1)*de+di
      ebtm=dble(e(1))
      return
      end
