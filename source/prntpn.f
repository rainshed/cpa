      subroutine prntpn(tchc,tchs,fcs,ew,ez,ng,mxl,meshe,isr)
c-----------------------------------------------------------------------
c     Given the Tchebycheff expansion tchc, and tchs and the
c     shoft parameter fcs, this program display s and c function.
c     coded by H.Akai, 1984, Juelich
c     adapted to spin-orbit version, 18 April, 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension tchc(ng,mxl**2),tchs(ng,mxl**2),tc(30),sn(5),cn(5)
     &         ,fcs(3,mxl**2)
      logical sra
      data c/274.0720442d0/
      sra=isr .eq. 1
      if(ng .gt. 30) call errtrp(1,'prntpn','ng too large')
      write(6,1100)
      d=2d0/dble(meshe-1)
      tc(1)=1d0
      do 20 k=1,meshe
      tc(2)=-1d0+d*dble(k-1)
      e=ez*tc(2)+ew
      if(e .lt. 0d0) go to 20
      erl=e
      if(sra) erl=e*(1d0+e/c**2)
      es=sqrt(erl)
      do 30 l=3,ng
   30 tc(l)=2d0*tc(2)*tc(l-1)-tc(l-2)
      do 40 j=1,mxl
      jj=j*(j-1)+1
      eff=e-fcs(3,jj)
      cn(j)=0d0
      sn(j)=0d0
      do 50 l=1,ng
      cn(j)=cn(j)+tchc(l,jj)*tc(l)
   50 sn(j)=sn(j)+tchs(l,jj)*tc(l)
      cn(j)=cn(j)*eff+fcs(1,jj)
      sn(j)=sn(j)*eff+fcs(2,jj)
      omeg=sqrt(cn(j)**2/erl**(j-1)+sn(j)**2*erl**j)
      cn(j)=cn(j)/es**(j-1)/omeg
   40 sn(j)=-sn(j)*es**j/omeg
      write(6,1200)k,e,(cn(j),j=1,mxl),(sn(j),j=1,mxl)
   20 continue
c     write(6,1300)(k,(tchc(k,l*(l-1)+1),l=1,mxl)
c    &             ,(tchs(k,l*(l-1)+1),l=1,mxl),k=1,ng)
c1300 format(/'   ===tchc,tchs ==='/(5x,i3,5x,1p,3e14.7,5x,3e14.7))
      return
 1000 format('   ***err in prpha...ng=',i3,' too large')
 1100 format(//t3,'k',t10,'x',t34,'c(x)',t69,'s(x)'/2x,75('-'))
 1200 format(1x,i3,f6.2,3x,3f10.5,3x,3f10.5)
      end
