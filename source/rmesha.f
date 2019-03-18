      subroutine rmesha(x1,x2,x3,dr,xr,meshr)
c-----------------------------------------------------------------------
c     +----------------------------------------------------+
c         radial mesh of the form r(x)=c*x**a*(x+g)**b+d
c     +----------------------------------------------------+
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension xr(meshr),dr(meshr)
      data a,b,g,e/6d0,-4d0,1d-1,-10d0/
      if(mod(meshr,2) .eq. 0) go to 20
      write(6,1000)meshr
 1000 format('   ***err in rmesha...odd meshr =',i5)
      stop
c
   20 r1=(1d0-e)/(dble(meshr-1)-e)
      c=(x1-x2)/(r1**a*(r1+g)**b-(1d0+g)**b)
      d=x2-c*(1d0+g)**b
      do 10 k=1,meshr-1
      x=(dble(k)-e)/(dble(meshr-1)-e)
      xr(k)=c*x**a*(x+g)**b+d
   10 dr(k)=(xr(k)-d)*(a/x+b/(x+g))/(dble(meshr-1)-e)
      xr(meshr)=x3
      dr(meshr)=x3-x2
c     write(6,1100)xr(1),xr(2),xr(meshr-2),xr(meshr-1)
c     write(6,1200)dr(1),dr(2),dr(meshr-2),dr(meshr-1)
c1100 format('   x1,x2,xmeshr-2,xmeshr-1=',1p,4e14.7)
c1200 format('   d1,d2,dmeshr-2,dmeshr-1=',1p,4e14.7)
      return
      end
