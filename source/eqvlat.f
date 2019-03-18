      function eqvlat(r,n,g,ftype,gsf,ng)
c----------------------------------------------------------------------
c     This function checks if the set of the n lattice points r in the
c     unit cell is totally equivalent to the set that is characterized
c     by its geometrical structure factors gms(g)=abs(sum{exp(igr)}),
c     that is, equivalent within the translation by the basic lattice
c     vector and any simultaneou shifts. Here g's are the ng different
c     reciprocal lattice vectors. Actually g's and r's are given in the
c     unit of 2pi/a and the lattice constant a, respectively.
c     Coded by H. Akai, 30 Nov 1996, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 r(3,n),g(3,ng),ftype(n),gsf(ng)
      complex*16 cunit,c
      logical eqvlat
      data zero/1d-5/
c     write(*,'((1x,i3,f12.5,2x,3f12.5))')
c    &    (j,ftype(j),(r(i,j),i=1,3),j=1,n)
      eqvlat=.false.
      twopi=8d0*atan(1d0)
      cunit=dcmplx(0d0,twopi)
c     --- loop over g's
      do 10 i=1,ng
      c=(0d0,0d0)
c     --- sum of exp(igr) over r's for a given g
      do 20 j=1,n
   20 c=c+ftype(j)*exp(cunit*(g(1,i)*r(1,j)
     &        +g(2,i)*r(2,j)+g(3,i)*r(3,j)))
c     write(11,'(1x,i3,2f10.7)')i,c
      if(abs(abs(c)-gsf(i))/dble(n) .gt. zero) return
   10 continue
      eqvlat=.true.
      end
