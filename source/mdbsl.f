      subroutine mdbsl(e,l,rj,dj,rn,dn)
c-----------------------------------------------------------------------
c         Modified Bessel & Neuman functions and their derivatives
c
c         j(z,l) = z**l * rj(l+1)
c         n(z,l) = z**(-l-1) * rn(l+1)    ( e = z**2, l=0,1,2... )
c
c         Definition of j and n is after Messiah, i.e. the sign of
c         the Neuman function is different.
c         coded by H.Akai, 1975, Osaka
c         revised by H.Akai, 1989, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (lmx=10, lmx2=lmx+2, lmx3=2*lmx+3)
      real*8 rj(*),dj(*),rn(*),dn(*),f(lmx3),g(lmx3),cj(lmx2)
      logical first
      save first,cj
      data first/.true./, small/1d-16/
      if(first) then
      first=.false.
      cj(1)=1d0
      do 10 j=2,lmx2
   10 cj(j)=cj(j-1)/dble(2*j-1)
      endif
      if(l .lt. 0 .or. l .gt. lmx) then
      write(6,1000)
 1000 format ('   ***err in mdbsl...illegal l')
      stop
      endif
      do 30 i=1,2
      j=l+2-i
      jj=2*j+1
      b=cj(j+1)
      a=b
      do 20 m=2,2000,2
      b=-b*e/dble(m*(m+jj))
      a=a+b
      if(abs(b) .lt. abs(a)*small) go to 30
   20 continue
   30 f(i)=a
      g(2)=dble(l)*f(2)-e*f(1)
      imx=2*l+1
      do 40 i=1,imx
      n=l+1-i
      f(i+2)=dble(2*n+1)*f(i+1)-e*f(i)
   40 g(i+2)=dble(n-1)*f(i+2)-e*f(i+1)
      do 50 i=1,l+1
      j=l+3-i
      rj(i)=f(j)
      dj(i)=g(j)
      j=l+2+i
      sgn=(-1d0)**(i-1)
      rn(i)=f(j)*sgn
   50 dn(i)=g(j)*sgn
      return
      end
