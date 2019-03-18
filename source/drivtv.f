      subroutine drivtv(f,kmx,xr,dr,k,df,df2)
c-----------------------------------------------------------------------
c     Given an array f(kmx), this program returns the first and second
c     derivatives 'df' and 'df2' of 'f', at a specified grid point 'k'.
c     The five point fomula is used for the 1st and 2nd erivative.
c     Coded by H. Akai 27 Nov. 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(morder=3,ndata=15)
      real*8 f(kmx),xr(kmx),dr(kmx),a(5,2),p(morder)
      data a/1d0,-8d0,0d0,8d0,-1d0
     &     ,-1d0,16d0,-30d0,16d0,-1d0/
c     if(kmx .lt. 5) call errtrp(1,'drvtv2','kmx < 5 is illegal')
      df=0d0
      df2=0d0
      ddr=0d0
      if(k .gt. 2 .and. k .lt. kmx-1) then
      do 10 j=1,5
      ddr=ddr+a(j,1)*dr(k-3+j)
      df=df+a(j,1)*f(k-3+j)
   10 df2=df2+a(j,2)*f(k-3+j)
      ddr=ddr/12d0
      df=df/12d0
      df2=df2/12d0
      df2=(-df*ddr/dr(k)+df2)/dr(k)**2
      df=df/dr(k)
      else
      if(k .le. 2) then
      call lsqftr(xr(1),f(1),ndata,p,morder)
      else
      call lsqftr(xr(kmx+1-ndata),f(kmx+1-ndata),ndata,p,morder)
      endif
      do 20 j=1,morder-2
      q=p(j)*dble(morder-j)
      df=df*xr(k)+q
   20 df2=df2*xr(k)+q*dble(morder-1-j)
      df=df*xr(k)+p(morder-1)
      endif
c     write(*,'(1x,1p,4e15.7)')xr(k),f(k),df,df2
      end
