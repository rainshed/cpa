      subroutine stfact(r,ftype,natm,g,gsf,ng)
c----------------------------------------------------------------------
c     This program calculates the geometrical structure factor gsf(g)
c     for given set of n atomic positions r. gsf(g) is given by
c     gsf(g)=abs{sum{ftype*exp(igr)}). Here g's are the 3n different
c     reciprocal lattice vectors. g's and r's are given in
c     the unit of 2pi/a and the lattice constant a, respectively.
c     Coded by H. Akai, 30 Nov 1996, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 r(3,natm),g(3,ng),ftype(natm),gsf(ng)
      complex*16 cunit,c
      twopi=8d0*atan(1d0)
      cunit=dcmplx(0d0,twopi)
c     write(*,*)'ng=',ng
c     --- loop over g's
      do 10 i=1,ng
      c=(0d0,0d0)
c     --- sum of exp(igr) over r's for a given g
      do 20 j=1,natm
   20 c=c+ftype(j)*exp(cunit*(g(1,i)*r(1,j)
     &          +g(2,i)*r(2,j)+g(3,i)*r(3,j)))
      gsf(i)=abs(c)
c     write(*,'(1x,i3,2f10.7)')i,c
   10 continue
c     write(11,'((1x,i3,3f12.5,2x,f12.5))')
c    &  (i,(g(j,i),j=1,3),gsf(i),i=1,ng)
      end
