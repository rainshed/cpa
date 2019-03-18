      subroutine dplini(v,corlvl,rstr,tm,ng,dr,xr,meshr,mmxl,wk,isr       
     &                 ,match,dpl,asa)
c------------------------------------------------------------------------
c     This program calculated square matrix elements |<1s| r |p(e)>|**2,
c     which are used in the calculation of the K-edge dipole transition
c     probability.
c     Coded by H. Akai, 20 May 1997, Osaka
c     Modified by H. Akai, Tokyo, 24 Jan. 2018.
c------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(meshr),corlvl(18),rstr(meshr,mmxl,ng),dr(meshr)
     &      ,xr(meshr),dpl(ng,2),tm(ng,ng),wk(meshr,2)
      integer match(18),indx(2)
      logical asa
c     --- indx points (l=1,m=-1) and (l=1,m=1) states.
      data indx/2,4/
      call clrarr(dpl,2*ng)
c     --- get 1s core wavefunction.
      call corada(corlvl(1),1,wk(1,1),rin,match(1),g1,g2,node,v,dr,xr
     &           ,meshr,isr,asa)
c     --- corada returns wk=R(r)**2, where R(r)=phi(r)*r is the radial
c         wavefunction.
      do 10 k=1,meshr-1
c     --- The 1s wavefunction has no node.
   10 wk(k,1)=sqrt(dabs(wk(k,1)))/xr(k)
      do 20 i=1,2
      ii=indx(i)
      do 20 j=1,ng
      do 30 k=1,meshr-1
c     --- rstr stores R(r)=phi(r)*r. In this case, the dipole matrix element
c         <1s|r|p(e)> is given by the following expression. 
   30 wk(k,2)=wk(k,1)*rstr(k,ii,j)
      s=fintgr(wk(1,2),dr,xr,meshr)
c     --- transion probability is propotional to s**2
      s=s**2
c     write(*,'(1x,a,f12.7)')'   matrix element=',s
      do 20 l=1,ng
   20 dpl(l,i)=dpl(l,i)+s*tm(l,j)
      end
      subroutine dplmtx(e,ew,ez,dipole,dpl,tc,ng)
c------------------------------------------------------------------------
c     Given an energy, this part calculates dipole matrix elements by
c     use of its Tchebycheff expansion coefficients obtained by the
c     first part of this program.
c------------------------------------------------------------------------
      real*8 tc(ng),dpl(ng,2),dipole(2)
      call gntcs(e,ew,ez,tc,ng)
      do 40 i=1,2
      dipole(i)=0d0
      do 40 j=1,ng
   40 dipole(i)=dipole(i)+tc(j)*dpl(j,i)
c     write(*,'(1x,a,3f12.7)')'e,dipole=',e,dipole
      end
