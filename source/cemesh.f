      subroutine cemesh(ef,ewidth,edelt,ebtm,e,kmx)
c--------------------------------------------------------------------
c     Generate a semielliptic energy contour. Mesh points are located
c     following fermi's distribution function such that they are
c     distributed densely near the real axis.
c     The following are the examples of typical cases used in the past.
c       edelt/3d-4/, r/0.65d0 /, h/2d-1/, rc/4d-1/
c       edelt/1d-3/, r/0.70d0 /, h/2d-1/, rc/4d-1/
c       edelt/1d-3/, r/1.20d0 /, h/2d-1/, rc/4d-1/
c       edelt/3d-3/, r/1.40d0 /, h/2d-1/, rc/4d-1/
c       edelt/1d-4/, r/0.40d0 /, h/5d-1/, rc/4d-1/
c     coded by H.Akai, 1983, Juelich
c     latest version, 30 Nov.1997, Osaka
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 e(kmx)
c     data h/5d-1/,rc/4d-1/,one/0.99999999d0/
      data h/2d-1/,rc/4d-1/,one/0.99999999d0/
      r=ewidth/2d0
      pi=4d0*atan(1d0)
      kc=rc*dble(kmx)+5d-1
      beta=log(pi/(edelt/h/r)-1d0)/dble(kmx-kc)
      f=pi*(exp(beta*dble(1-kc))+1d0)*one
      do 10 k=1,kmx
      theta=f/(exp(beta*dble(k-kc))+1d0)
   10 e(k)=ef+r*dcmplx(cos(theta)-1d0,h*sin(theta))
      ebtm=dble(e(1))
      return
      end
