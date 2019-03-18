      subroutine tchstr(eta,vc,nkp,vectk,lmx,np,nt,vset,tch
     &                 ,e1,e2,ng,mch,rpt,dr,nrpt,pexf)
c-----------------------------------------------------------------------
c     Tchebycheff expansion of the structural Green function.
c     In order to avoid overflows that often occurs in the
c     series expansion, now the coeeficients gieven by
c     "mseque" involves the factorial factors and energy
c     normalization factors. This version is modified so as
c     adapted to this modification (7 July 2007.)
c     Coded by M. and H.Akai, 1980, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c     Revised by H. Akai, 7 July 2007, Osaka
c     Modified by H. Akai to adapt OpenMP version, 6 Jan 2016, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     --- in the parameter, lx=3 means that the highest anguler states
c         that are considered are f type, and so on.
      parameter (nd=300)
      complex*16 tch(mch,ng,nkp),pexf(np,nkp),cexf
      real*8 rpt(3,nrpt),dr(3),vectk(3,nkp),vset(3,nt)
     &         ,vectr(3),vect(3)
      complex*16,allocatable::gexf(:),rexf(:),d(:),phf(:)
      real*8,allocatable::tc(:,:),s(:,:,:),sy(:,:),d3(:),h(:,:)
     &      ,pr(:),ak(:)
      integer,allocatable::lofj(:)
      data small/1d-3/, expmin/0d0/
      ntch=(2*lmx+1)**2
      mx2=2*lmx+1
      allocate(gexf(nt),rexf(nrpt),d(ntch),phf(mx2),tc(ng,ng)
     &        ,s(mx2,ng,nrpt),sy(ntch,nrpt),d3(ng),h(ntch,nt)
     &        ,pr(nt),ak((nd+1)*mx2),lofj(ntch),stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'tchstr','allocation fails')
      pi=4d0*atan(1d0)
      jmax=lmx
      jjmax=(jmax+1)**2
c     --- some constants
c     --- oldeta is the eta used in the old version.
      twoeta=2d0*eta
      oldeta=twoeta**2
      o1=-4d0*eta
      q1=2d0/pi/vc
      q2=-4d0*sqrt(pi)
      p1=pi/(2d0*dble(ng))
      p2=2d0/dble(ng)
      ew=(e2+e1)*5d-1
      ez=(e2-e1)*5d-1
      en=max(abs(expmin),abs(e2))
      en=max(en,1d0)
c     --- phase factor
      do 10 j=1,jmax+1
   10 phf(j)=(0d0,-1d0)**(j-1)
c     --- mapping (l,m) --> l
      do 20 j=1,jjmax
   20 lofj(j)=sqrt(dble(j-1))+1.000001d0
c     --- generate Tchebycheff polinomials
      do 30 k=1,ng
      ang=2d0*dble(k-1)+1d0
      ang=p1*ang
      tc(1,k)=1d0
      tc(2,k)=cos(ang)
      do 30 l=3,ng
   30 tc(l,k)=2d0*tc(2,k)*tc(l-1,k)-tc(l-2,k)
c
c     --- first prepair quantities which do not depend on k.
c
c     --- d3 depends only on energy
      if(dr(1)**2+dr(2)**2+dr(3)**2 .lt. small) then
c     --- loop over energies
      do 40 k=1,ng
      e=ew+ez*tc(2,k)
      o2=e/oldeta
c     --- fd3 returns d3(e/eta)/o1
   40 d3(k)=fd3(o2)*o1
c     --- in the lattice diagonal case (0,0,0) terms in the lattice
c     --- sum appearing below should be omitted.
      n1=2
      else
C     --- d3 disappears in the site off diagonal cases.
      do 50 k=1,ng
   50 d3(k)=0d0
c     --- but (0,0,0) should be included in the lattice sum.
      n1=1
      endif
c
c     --- gaussian integral for all l's, e's and the lattice points.
c     --- loop over r
      do 60 n=n1,nrpt
      do 70 i=1,3
   70 vectr(i)=pi*(rpt(i,n)-dr(i))
      da=sqrt(vectr(1)**2+vectr(2)**2+vectr(3)**2)
      call mseque(jmax,ak,da,eta,nd,en)
      call realh(vectr,sy(1,n),jmax)
c     --- loop over energies
      do 60 k=1,ng
      e=ew+ez*tc(2,k)
      if(e .gt.  expmin) then
      enrm=e/en
c
c     --- series expansion in the form of a recursive relation
c     --- is used for  e>-2. A direct integratin becomes very
c     --- hard for e>0.
      do 80 m=1,jmax+1
      m0=(nd+1)*(m-1)
      s(m,k,n)=0d0
      do 80 i=nd+1,1,-1
   80 s(m,k,n)=s(m,k,n)*enrm+ak(i+m0)
      else
c
c     --- numerical integration for e<0 cases. This treatment is
c     --- needed since the above series hardly gives an accurate
c     --- result for e<0 due to round-off error (-2<e<0 is still
c     --- safe if 'erfc' is accurate enough), while the
c     --- integration becomes rather easy in this case.
      do 90 m=1,jmax+1
      mone=m-1
c     write(*,*)'n,k,m=',n,k,m
      call qromo2(da,e,mone,twoeta,1d60,s(m,k,n))
   90 continue
      endif
   60 continue
c
      call clrarc(tch,mch*ng*nkp)
c     --- big loop over k points
      do 110 kp=1,nkp
c     --- g loop
      do 120 i=1,nt
      do 130 j=1,3
  130 vect(j)=vectk(j,kp)+vset(j,i)
      call realh(vect,h(1,i),jmax)
      xx=2d0*pi*(vect(1)*dr(1)+vect(2)*dr(2)+vect(3)*dr(3))
      gexf(i)=dcmplx(cos(xx),sin(xx))
  120 pr(i)=vect(1)**2+vect(2)**2+vect(3)**2
c     --- store (k+g)**2, yl and exp(i*(k+g)*dr) for preferred set
      do 140 i=1,np
  140 pexf(i,kp)=gexf(i)
c     --- r loop
      do 150 i=n1,nrpt
      xx=2d0*pi*(vectk(1,kp)*rpt(1,i)+vectk(2,kp)*rpt(2,i)
     &  +vectk(3,kp)*rpt(3,i))
  150 rexf(i)=q2*dcmplx(cos(xx),sin(xx))
c     --- energy loop
      do 160 k=1,ng
      e=ew+ez*tc(2,k)
      d(1)=dcmplx(d3(k),0d0)
      do 170 i=2,jjmax
  170 d(i)=(0d0,0d0)
      do 180 j=1,jjmax
      l=lofj(j)
      do 180 i=n1,nrpt
  180 d(j)=d(j)+rexf(i)*sy(j,i)*s(l,k,i)*phf(l)
c     --- loop over preferred set
      do 190 i=1,np
      x=(e-pr(i))/oldeta
      if(abs(x) .lt. small) then
      cc=q1*(1d0+(x/2d0)*(1d0+(x/3d0)*(1d0+(x/4d0))))/oldeta
      else
      cc=q1*(exp(x)-1d0)/(e-pr(i))
      endif
      cexf=cc*gexf(i)
      do 190 ml=1,jjmax
  190 d(ml)=d(ml)+cexf*h(ml,i)
c     --- loop over remainder set
      do 200 i=np+1,nt
      cc=q1*exp((e-pr(i))/oldeta)/(e-pr(i))
      cexf=cc*gexf(i)
      do 200 ml=1,jjmax
  200 d(ml)=d(ml)+cexf*h(ml,i)
      do 160 ml=1,jjmax
      do 160 l=1,ng
  160 tch(ml,l,kp)=tch(ml,l,kp)+d(ml)*tc(l,k)
      do 110 ml=1,jjmax
      tch(ml,1,kp)=5d-1*p2*tch(ml,1,kp)
      do 110 l=2,ng
  110 tch(ml,l,kp)=p2*tch(ml,l,kp)
      deallocate(gexf,rexf,d,phf,tc,s,sy,d3,h,pr,ak,lofj)
      return
      end
