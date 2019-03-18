      subroutine srtrns(u,mxl)
c---------------------------------------------------------------------
c     Generate the unitary transformation matrix u that transforms
c     spherical harmonics to real harmonics. mxl=l+1 where l is the
c     maximum angular momentum. The transformation must be performed
c     wtihin each subspace of l composed of 2*l+1 bases.
c     real Y = conjg(u) * spherical Y
c     When u is used for matrix transformation from spherical to
c     real harmonics representation,
c     (real representaion) = u * (sperical representaion) * u^+
c     must be used.
c     In this program, u=u(i,j,1) and u^+=conjg(u(j,i,2)=u(j,i,2)
c     Coded by H. Akai, 9 April 2014, Tokyo
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 u(2*mxl-1,2*mxl-1,2),cunit
      data cunit/(0d0,1d0)/
      a=sqrt(2d0)/2d0
      do 10 j=1,2*mxl-1
      do 10 i=1,2*mxl-1
   10 u(i,j,1)=(0d0,0d0)
c     --- construct u, where  (real y) = conjg(u) * (spherical y)
      u(mxl,mxl,1)=(1d0,0d0)
      do 20 m=1,mxl-1
      sgn=(-1d0)**m
      mr=mxl-m
      mc=mxl+m
      u(mr,mr,1)=-a*cunit
      u(mr,mc,1)=sgn*a*cunit
      u(mc,mr,1)=a
   20 u(mc,mc,1)=sgn*a
      do 30 j=1,2*mxl-1
      do 30 i=1,2*mxl-1
   30 u(i,j,2)=conjg(u(j,i,1))
      end
