      subroutine mseque(lmax,ak,da,eta,nd,en)
c-----------------------------------------------------------------------
c     Calculate integ( x**(2*l)*exp(-(x*r)**2) dx )
c     da=pi*r
c     en is an energy normalization factor.
c     See H.A. note.
c
c     In order to avoid overflows the factorials are
c     included in the coefficietns. Also the energies
c     are assumed to be normalized by the maximam value
c     of the energy.
c     Coded by M. and H.Akai, 1980, Osaka
c     Major revision by H. Akai on 7 July 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (small=1d-70)
      real*8 ak(0:nd,0:lmax)
      pi=4d0*atan(1d0)
      x=2d0*eta*da
      fc=2d0*da**2
      ak(0,0)=sqrt(pi)*erfc(x)/(2d0*da)
      ext=exp(-x**2)
      do 10 l=1,lmax
   10 ak(0,l)=(dble(2*l-1)*ak(0,l-1)+(2d0*eta)**(2*l-1)*ext)/fc
      do 20 l=0,lmax
      p=1d0
      do 20 i=1,nd
      if(abs(ak(i-1,l)) .gt. small) then
      p=p*en/dble(i)
      if(abs(p) .lt. small) p=0d0
      ak(i,l)=(en*fc*ak(i-1,l)/dble(i)-p*(2d0*eta)**(2*(l-i)+1)*ext)
     &      /dble(2*(l-i)+1)
      else
      ak(i,l)=0d0
      endif
   20 continue
c     write(*,*)(ak(i,0),i=0,nd)
      return
      end
