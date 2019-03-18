      subroutine cpshft(e,tchc,tchs,cm,fcs,gs,tc,mxl,mxlcmp,ng,t,pf,str
     &                 ,isr,ids)
c-----------------------------------------------------------------------
c     In this version, pf does not include the phase factor i^(l1-l2)
c     any more. It is now attached to the scattering path operator
c     calculated  by kkrsed.
c     Modified by H. Akai, 30 July 2005, Osaka
c
c     ------------------------------------
c     --- KKR-CPA + spin-orbit version ---
c     ------------------------------------
c     Given Tchebycheff expansion coeficients tchc,tchs and the
c     shift constants fcs, this program returns t-matrix (or its
c     inverse) and the phase factor fanction, t and pf.
c     To facilitate the calculation of the phase of the Green's
c     function this program also returns function 'str'.
c     coded by H.Akai, 1983, Juelich
c     Only the modification is that now 'str' is defined differently.
c     KKR CPA implemented by H.Akai, 21 Sep 1996, Osaka
c     adapted to spin-orbit version, 20 April, 1997, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 tc(ng),t(mxl**2),pf(mxl**2,mxl**2),str(mxl**2)
     &          ,p(36),gs(mxl**2),e,ek,ekl,ea,erel,eef,ex,cx,sx,fx
      real*8 tchc(ng,mxl**2),tchs(ng,mxl**2),fcs(3,mxl**2)
     &      ,cm(ng,mxl**2)
      data c/274.0720442d0/
      mmxl=mxlcmp**2
      erel=e
      if(isr .eq. 1) erel=e+(e/c)**2
      ek=sqrt(erel)
      do 10 l=1,mmxl
      ll=sqrt(dble(l)-5d-1)+1
      eef=e-fcs(3,l)
      fx=(0d0,0d0)
      cx=(0d0,0d0)
      sx=(0d0,0d0)
      do 20 n=1,ng
      fx=fx+cm(n,l)*tc(n)
      cx=cx+tchc(n,l)*tc(n)
   20 sx=sx+tchs(n,l)*tc(n)
      cx=cx*eef+fcs(1,l)
      sx=sx*eef+fcs(2,l)
      ekl=ek**(ll-1)
      ea=erel**(ll-1)
      ex=cx+(0d0,1d0)*sx*ea*ek
c     write(*,'(1x,i2,7f10.5)')l,dble(e),ex,cx-(0d0,1d0)*sx*ea*ek
      p(l)=ekl/ex
      str(l)=log(sx*ea)
c     str(l)=log(ex)
      gs(l)=(fx-(0d0,1d0)*ek*ekl*p(l))/cx
   10 t(l)=-cx/(sx*ea)-(0d0,1d0)*ek
c     --- check for spin-orbit coupling
c     do 40 l=1,mmxl
c     ll=sqrt(dble(l)-5d-1)+1
c     lc=(ll-1)**2+ll
c     t(l)=t(lc)
c     p(l)=p(lc)
c  40 gs(l)=gs(lc)
c     ---in this case, t(l) is the inverse t-matrix.
c  10 t(l)=-sx*p(l)*ekl
c     ---in this case, t(l) is the t-matrix.
      if(ids .eq. 2) write(*,'(1x,7f10.5)')dble(e),gs(1),gs(3),gs(7)
c     write(*,'(1x,5f10.5)')dble(e),cx,sx*ea*ek
c     write(*,'(1x,7f10.5)')dble(e),t(1),t(2),t(3)
c     ---gs(l) is the single potential Green's function.
c     ---str(l) is the quantity needed to calculate
c     ---phase of the KKR Green's function.
c        The phase factor i**(l2-l3) is needed since it is
c        not attached in the subroutine kkrsed. In the
c        definitioan of the Green function it should be
c        attached though it is rather harmless untill
c        the very end where the phase of the wave fucntion
c        is considered.
      do 30 l1=1,mmxl
      ll1=sqrt(dble(l1)-5d-1)+1
      do 30 l2=1,mmxl
      ll2=sqrt(dble(l2)-5d-1)+1
   30 pf(l1,l2)=p(l1)*p(l2)*(1d0,0d0)**(ll1-ll2)
      return
      end
