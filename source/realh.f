      subroutine realh(v,y,lmax)
c-----------------------------------------------------------------------
c     Given a vector v(3) and lmax this program calculoates r**l times
c     real spherical hermonics, in the order of y(0), y(1,sin-type),
c     y(1,0),y(1,cos-type),...,y(lmax,sin-type),..,y(lmax,0),
c     y(lmax,cos-type).
c     coded by M.Akai and H.Akai, 1980, Osaka
c     latest revision on April 1, 1992
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (lmx=12,lp=lmx+1,lb=2*lmx+1)
      dimension v(3),y(*),a(lb),b(lp,lp),t(lp,lp),pog(lp),c(lp),s(lp)
      logical init
      save b,pog,init
      data init/.true./
c
c   data for b and pog are generated once and all at the first call.
c   lmx in parameter statement should be changed for larger lmax.
      if(init) then
      init=.false.
      pi=4d0*atan(1d0)
      a(1)=1d0
      dd=1d0
      do 10 i=2,lb
      dd=dd*dble(i-1)
   10 a(i)=sqrt(dd)
      b(1,1)=5d-1/sqrt(pi)
      do 20 l=2,lp
      x=dble(2*l-1)/pi
      b(l,1)=sqrt(x)/2d0
      do 20 n=2,l
      k1=l-n+1
      k2=l+n-1
   20 b(l,n)=sqrt(x/2d0)*a(k1)/a(k2)
      pog(1)=1d0
      do 30 l=2,lp
   30 pog(l)=pog(l-1)*dble(2*l-3)
      endif
c
      if(lmax .lt. 0 .or. lmax .gt. lmx) then
      write(6,1000)lmax
 1000 format('   ***err in realh...lmax=',i4,' is illegal')
      stop
      endif
c
      x3=v(1)**2+v(2)**2
      r=x3+v(3)**2
      la=lmax+1
      do 40 i=1,la
   40 t(i,i)=pog(i)
      c(1)=1d0
      s(1)=0d0
      do 50 n=2,la
      c(n)=v(1)*c(n-1)-v(2)*s(n-1)
   50 s(n)=v(1)*s(n-1)+v(2)*c(n-1)
      t(2,1)=v(3)
      do 60 l=3,la
   60 t(l,1)=(dble(2*l-3)*v(3)*t(l-1,1)-r*dble(l-2)*t(l-2,1))/dble(l-1)
      if(v(1)**2+v(2)**2 .lt. v(3)**2) then
c---  following unneccesary substitution is to fool the nec compiler
      v3=v(3)
      do 70 l=3,la
      do 70 n=2,l-1
      nn=l+1-n
   70 t(l,nn)=(r*dble(l+nn-2)*t(l-1,nn)-x3*t(l,nn+1))/dble(l-nn)/v3
      else
      do 80 l=3,la
      do 80 n=2,l-1
   80 t(l,n)=(r*dble(l+n-3)*t(l-1,n-1)-v(3)*dble(l-n+1)*t(l,n-1))/x3
      endif
      do 90 l=1,la
      do 90 n=1,l
   90 t(l,n)=b(l,n)*t(l,n)
      y(1)=t(1,1)
      i=1
      if(lmax .eq. 0) return
      do 110 l=2,la
      do 100 n=2,l
      nn=l-n+2
      i=i+1
  100 y(i)=t(l,nn)*s(nn)
      i=i+1
      y(i)=t(l,1)
      do 110 n=2,l
      i=i+1
  110 y(i)=t(l,n)*c(n)
      return
      end
