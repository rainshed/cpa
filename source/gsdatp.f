      subroutine gsdatp(wk,a,vc,atmicp,r,anclr,corl,ro,rmt,dr,xr
     &                 ,iwk,itype,natm,ntyp,meshr,ncmp,ncmpx,conc)
c-----------------------------------------------------------------------
c     -----------------------
c     --- KKR-CPA version ---
c     -----------------------
c     Initial potential for self-consistent calculation of
c     electronic structures of metals and compounds with
c     a complex lattice.
c     Mattheiss' prescription is used. Atomic potentials
c     are generated self-consistently by LSD scheme.
c     No spin density is included.
c
c     Work area of size 8*mxa+nb of real*8 and 2*nb of integer needed.
c
c     coded by H.Akai on Jan. 16, 1986 (Osaka)
c     revised by H.Akai April, 1992, Osaka
c     KKR-CPA implemented by H.Akai, 7 Sep. 1996, Osaka
c     Modified by H. Akai, Tokyo, Jan, 26, 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(nb=1200,mxa=441)
c     mxa is radial mesh number used in atomic potential program
c     and nb is the maxumum number of neighboring atoms to be
c     considered.
      real*8 wk(8*mxa+2*nb),atmicp(3,natm),r(3,3),cnf(18)
     &      ,anclr(ncmpx),corl(18,ncmpx),conc(ncmpx)
     &      ,ro(meshr,ncmpx),rmt(ntyp)
     &      ,dr(meshr,ncmpx),xr(meshr,ncmpx)
      integer iwk(16*mxa+6*nb),itype(natm),ncmp(ntyp)
      integer,allocatable::iwt(:)
      logical*1,allocatable:: do(:)
      data zero/1d-10/ , small/1d-7/
      allocate(iwt(ncmpx),do(ncmpx))
c     ---pointer ipds, ipro, ipda, ipxa, iptyp, etc. to the work
c        area are used. each corresponds to dist, ro, da, xa, ityp,
c        and so on. keep enough space for these regions.
c
      pi=4d0*atan(1d0)
c
      ipro=1
      ipda=mxa+1
      ipxa=2*mxa+1
      ipz0=3*mxa+1
      ipz1=4*mxa+1
      ipz2=5*mxa+1
      ipwk=6*mxa+1
      iprc=7*mxa+1
      ipdist=8*mxa+1
      ipwt=8*mxa+nb+1
      iptyp=2*(8*mxa+2*nb)+1
      ipdeg=2*(8*mxa+2*nb)+nb+1
c     iptyp=1
c     ipdeg=nb+1
c
c     ---'ratm' corresponds to the average atomic volume.
      atvol=a**3*vc/dble(natm)
      ratm=(3d0*atvol/4d0/pi)**(1d0/3d0)
c     ---define a radial mesh for each site and each component.
      do 10 i=1,ntyp
      do 10 j=1,ncmp(i)
      call jip(i,j,ji)
      iwt(ji)=0
      do(ji)=.true.
      rtin=a*rmt(i)
   10 call rmesha(1d-6,rtin,ratm,dr(1,ji),xr(1,ji),meshr)
c     ---how often each type appears?
      do 20 i=1,natm
   20 iwt(itype(i))=iwt(itype(i))+1
      call clrarr(ro,meshr*ncmpx)
c     ---now construct charge densities.
      do 30 i=1,ntyp
      do 30 ic=1,ncmp(i)
      call jip(i,ic,ici)
      if(do(ici)) then
      call clrarr(cnf,18)
c     ---first perform atomic calculation.
      call atmicv(anclr(ici),cnf,corl(1,ici),wk(ipro),wk(ipda)
     &           ,wk(ipxa),wk(ipz0),wk(ipz1),wk(ipz2),wk(ipwk)
     &           ,wk(iprc),mxa,ier)
      do 60 k=mxa,2,-1
      kdlim=k
      if(wk(ipro+k-1) .gt. small) go to 70
   60 continue
   70 dlim=wk(ipxa+kdlim-1)/a
c     ---then perform matthiss' prescription, i.e., accumulate
c        all the charges extended over the different sites.
      do 40 j=i,ntyp
      id=j
      do 40 jc=1,ncmp(j)
      call jip(j,jc,jcj)
      if(abs(anclr(jcj)-anclr(ici)) .lt. zero) then
      do(jcj)=.false.
      jcmp=jc
      do 80 l=1,18
   80 corl(l,jcj)=corl(l,ici)
c     ---make map of near neighbors.
c     Unfortunately, the following routine mkemap may be unnecessarily
c     called more than once only to bring about the same data. 
c     Presently, I have not devised any clever way to avoid this.
      call mkemap(id,atmicp,itype,r,natm,dlim,wk(ipdist)
     &           ,iwk(iptyp),iwk(ipdeg),nb,nd)
      do 50 n=1,nd
      l=iwk(iptyp+n-1)
   50 wk(ipwt+n-1)=conc(jcj)*dble(iwk(ipdeg+n-1)*iwt(j))/dble(iwt(l))
      call cryrho(wk(ipro),wk(ipda),wk(ipxa),mxa,a,ro,xr,meshr,nd
     &           ,wk(ipdist),iwk(iptyp),wk(ipwt),wk(ipwk)
     &           ,jcmp,ncmp,ntyp,ncmpx)
      endif
   40 continue
      endif
   30 continue
c     do 90 i=1,ncmpx
c     s1=fintgr(ro(1,i),dr(1,i),xr(1,i),meshr)
c  90 write(*,'(1x,a,i3,2f12.5)')
c    &  'component,Z,charge=',i,anclr(i),s1
      end
