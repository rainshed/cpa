      subroutine prntcs(tchc,tchs,fcs,ew,ez,ng,mxl,meshe,isr)
c-----------------------------------------------------------------------
c     Given the Tchebycheff expansion tchc, and tchs and the
c     shoft parameter fcs, this program display s and c function.
c     Also the phase shift can be displayed.
c     coded by H.Akai, 1984, Juelich
c     adapted to spin-orbit version, 18 April, 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 tchc(ng,mxl**2),tchs(ng,mxl**2),fcs(3,mxl**2),ph(5)
      complex*16 e,erel,es,eef,tc(30),cn(5),sn(5),ex,omeg
      logical sra
      data c/274.0720442d0/,del/1d-3/
      pi=4d0*atan(1d0)
      sra=isr .eq. 1
      if(ng .gt. 30) call errtrp(1,'prntcs','ng too large')
      write(6,1100)
      de=2d0*ez/dble(meshe-1)
      e0=ew-ez
      do 20 k=1,meshe
      e=e0+de*dble(k-1)+(0d0,1d0)*del
      call cgntcs(e,ew,ez,tc,ng)
      erel=e
      if(sra) erel=e+(e/c)**2
      es=sqrt(erel)
      do 40 j=1,mxl
c     jj=j*(j-1)+1
      jj=j**2
      eff=e-fcs(3,jj)
      cn(j)=(0d0,0d0)
      sn(j)=(0d0,0d0)
      do 50 l=1,ng
      cn(j)=cn(j)+tchc(l,jj)*tc(l)
   50 sn(j)=sn(j)+tchs(l,jj)*tc(l)
      cn(j)=cn(j)*eff+fcs(1,jj)
      sn(j)=sn(j)*eff+fcs(2,jj)
      omeg=sqrt(cn(j)**2/erel**(j-1)+sn(j)**2*erel**j)
      cn(j)=cn(j)/es**(j-1)/omeg
      sn(j)=-sn(j)*es**j/omeg
      ex=cn(j)+(0d0,1d0)*sn(j)
      ph(j)=dimag(log(ex))
      if(ph(j) .lt. -5d-1*pi) ph(j)=ph(j)+pi
      if(ph(j) .gt. 5d-1*pi) ph(j)=ph(j)-pi
   40 continue
      write(6,1200)k,dble(e),(dble(cn(j)),j=1,3),(dble(sn(j)),j=1,3)
c     write(6,1200)k,dble(e),(ph(j),j=1,mxl)
   20 continue
c     write(6,1300)(k,(tchc(k,l*(l-1)+1),l=1,mxl)
c                  ,(tchs(k,l*(l-1)+1),l=1,mxl),k=1,ng)
c1300 format(/'   ===tchc,tchs ==='/(5x,i3,5x,1p,3e14.7,5x,3e14.7))
      return
 1100 format(//t3,'k',t10,'e',t34,'c(e)',t69,'s(e)'/2x,75('-'))
 1200 format(1x,i3,f6.2,3x,3f10.5,3x,3f10.5)
c1100 format(//t3,'k',t8,'e',t20,'phase(e)'/2x,60('-'))
c1200 format(1x,i3,f6.2,3x,6f10.5)
      end
