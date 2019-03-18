      subroutine spmain
     & (wk,wkc,iwk,nwk,niwk,e,wt,detl,ff,tmt,phf,str,tc
     & ,anclr,rmt,ancp,field,atmicv,atmicp,elvl,tm,clks,zdmy
     & ,corlvl,dr,xr,amdlng,rms,v1,v2,v3
     & ,ro,rorg,cm,tchc,tchs,fcs,rstr,dosef,esic,sm,anc,tof
     & ,f,rmtd,q,hhf,match,itype,ncub,iblk,itblk
     & ,iwtyp,config,go,file,brvtyp,reltyp,sdftyp,magtyp,outtyp
     & ,type,atmtyp,cwtyp,title,bzqlty,record,ef0,dex,emxr,xlim
     & ,conc,mse,ng,mxl,tol,ids,inv,a,coa,boa,alpha,beta,gamma
     & ,edelt,ewidth,maxitr,pmxtyp,xmd,ntyp,natm,meshr,ndmx
     & ,lastmx,ncmp,ncmpx,gfree,ess
     & ,tcpa,phase,convrg,urotat,uu,irotat,ck,cj,iatm
     & ,lmxtyp,mxlcmp,msiz,lmxblk,iatmp,r,lmpair,lmpmx,angl,openc)
c-----------------------------------------------------------------------
c     ----------------------------------
c     --- KKR-CPA spin-orbit version ---
c     ----------------------------------
c     Actually this is the main program of kkr band structure
c     calculation. The program uses kkr-package developed since 1979.
c     Coded by H.Akai, 1986, Juelich
c     Latest revision, 3 Feb. 1996, Osaka
c     KKR-CPA implemented by H.Akai, 7 Sep. 1996, Osaka
c     Spin-Orbit version, 18 April, 1997, Osaka
c     Fixed spin moment procedure added, 8 Aug. 1999, Duisburg
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 e(mse,2),wt(ng,3,mse,2)
     &          ,wkc(mxl**2,ncmpx,mse),detl(mse,2),gfree(mse)
     &          ,ff(mse*mxl**4*natm)
     &          ,tmt(mse*mxl**2*ncmpx),phf(mse*mxl**4*ncmpx)
     &          ,str(mse*mxl**2*ncmpx),tc(mse*ng),zgiven
     &          ,ess(mse*(2*mxl-1)**2),tcpa(mse*mxl**4*natm,2)
     &          ,phase(mse*mxl**2*ncmpx)
     &          ,ck(mse*mxl**4*natm),cj(mse*mxl**4*natm)
     &          ,urotat((2*mxl-1)**2,mxl,24)
     &          ,uu(2*mxl-1,2*mxl-1,2)
      complex*16,allocatable::tch(:),pexf(:),spctrl(:)
c
      real*8     anclr(ncmpx),rmt(ntyp),ancp(ncmpx,2),field(ntyp)
     &          ,atmicp(3,natm),elvl(ng),tm(ng,ng),clks(lastmx)
     &          ,zdmy(ncmpx),corlvl(18,ncmpx,2),dr(meshr,ncmpx)
     &          ,xr(meshr,ncmpx),amdlng(natm,natm),rms(ncmpx,2)
     &          ,v1(meshr,ncmpx,2),v2(meshr,ncmpx,2),v3(meshr,ncmpx,2)
c
      real*8     wk(nwk),ro(meshr,ncmpx,2),dmy(50)
     &          ,rorg(20,ncmpx,2),tchc(ng,mxl**2,ncmpx)
     &          ,tchs(ng,mxl**2,ncmpx),fcs(3,mxl**2,ncmpx)
     &          ,rstr(meshr*mxl**2*ng,ncmpx,2)
     &          ,dosef(mxl**2,ncmpx,2),esic(ncmpx,2)
     &          ,sm(ng,mxl**2,ncmpx,2),anc(ncmpx)
     &          ,tof(mxl**2,ncmpx,2),f(mxl**2,ncmpx,ng,2)
     &          ,rmtd(natm),q(ncmpx),cm(ng,mxl**2,ncmpx)
     &          ,config(18,ncmpx,2),conc(ncmpx)
     &          ,hhf(ncmpx),xmd(ncmpx,mse,2,2)
     &          ,ebtm(2),r(3,3),g(3,3),ef(2),efs(2),tend(10),total(2)
     &          ,def(2),bnd2(2),protat(9,48),dipole(2),sftef(2),angl(3)
     &          ,gcornr(3,24),protdm(9,48,2),uprot(9,2,24)
      real*8,allocatable::vkp(:),wtkp(:),prr(:),rpt(:),gpt(:),hh(:)
     &          ,dpl(:,:),fm(:),wg(:),um(:),pm(:),gsf(:)
c
      integer    itype(natm),ncub(lastmx),iwk(niwk),ncmp(ntyp)
     &          ,iblk(natm,natm),iwtyp(ntyp)
     &          ,match(18,ncmpx,2),itblk(5,natm**2)
     &          ,irotat(natm*48),iatm(ntyp),lmxblk(ndmx)
     &          ,isymop(48),lmxtyp(ntyp),mxlcmp(ncmpx),iatmp(natm)
     &          ,lmpair(2,lmpmx),kcrt(25),isymdm(48,2)
      integer,allocatable::korder(:)
c
      character  type(ntyp)*8,atmtyp(natm)*8,atmicv(3,natm)*24
     &          ,go*3,file*256,brvtyp*6,reltyp*6,sdftyp*12,magtyp*4
     &          ,outtyp*6,title*400,token*80
     &          ,bzqlty*8,record*4,bravai*3,status*2
     &          ,cwtyp*4,today*11,logfil*265,inffil*265,spcfil(2)*265
     &          ,pmxtyp*16,trmkey*16
      logical    convrg(mse),cpacnv(2),msg,readin,alphac,asa,bckup
     &          ,ifkey,ereset,fbz,openc
      data       bckup/.true./
c
c
c     ------------------------------------------------------------------
c
      logfil=file//'.log'
      idmy=lftfil(logfil)
      inffil=file//'.info'
      idmy=lftfil(inffil)
      spcfil(1)=file//'_up.spc'
      spcfil(2)=file//'_dn.spc'
      idmy=lftfil(spcfil(1))
      idmy=lftfil(spcfil(2))
c     --- write on log-file ---
      if(go .eq. 'log') then
      call udate(today)
      open(26,file=logfil,form='formatted',
     &     status='unknown')
      do 97 irec=1,10000
   97 read(26,*,end=96)
      call errtrp(1,'spmain','too many records')
   96 backspace(26)
      write(26,'(1x,a)')today
      close(26)
      endif
c
      pi=4d0*atan(1d0)
      lmax=mxl-1
      mmxl=mxl**2
      pmix=rnmrdt(pmxtyp)
      write(*,1000)meshr,mse,ng,mxl
 1000 format('    meshr   mse    ng   mxl'/3x,4i6//)
c
c     ---print input data
      write(*,1900)go,file,brvtyp,a,coa,boa,alpha,beta,gamma
     &            ,edelt,ewidth,reltyp,sdftyp,magtyp,record,outtyp
     &            ,bzqlty,maxitr,pmix,ntyp,natm,ncmpx
 1900 format(/'   data read in'/
     &       '   go=',a3,'  file=',a64/'   brvtyp=',a6,
     &       '  a=',f9.5,'  c/a=',f8.5,'  b/a=',f8.5/
     &       '   alpha=',f5.1,'  beta=',f5.1,'  gamma=',f5.1/
     &       '   edelt=',1pe8.1,'  ewidth=',0pf7.3,'  reltyp=',a5,
     &       '  sdftyp=',a9/'   magtyp=',a4/'   record=',a4,
     &       '  outtyp=',a6,'  bzqlty=',a8/'   maxitr=',i3,
     &       '  pmix=',f8.5/'   ntyp=',i3,'  natm=',i3,
     &       '  ncmpx=',i3//)
c
      asa=ifkey('asa',sdftyp)
      if(asa) call getorg(eorg)
      call drvmsh(0d0,ewidth,edelt,ebtm(1),e(1,1),mse,ids)
      er=5d-1*ewidth
      eunder=-ebtm(1)
      write(*,2600)(k,e(k,1),k=1,mse)
 2600 format('   complex energy mesh'
     &      /(1x,3(i3,'(',f8.4,',',f7.4,') ')))
      write(*,'(1x)')
c
c     ---prepair i/o file
      call getfil(24,file)
c
      ls=0
      if(ifkey('ls',reltyp)) ls=1
      if(ifkey('nrl',reltyp)) then
      isr=0
      elseif(ifkey('sra',reltyp)) then
      isr=1
      else 
      call errtrp(1,'spmain','illegal reltyp')
      endif
      ibrav=ibrava(trmkey('tlt',brvtyp))
c     call ty2ity(atmtyp,natm,type,ntyp,itype)
c     if(asa) then
c     call chkasa(ibrav,a,coa,boa,alpha,beta,gamma,vc,r,g
c    &        ,atmicv,atmicp,itype,natm,anclr,rmt,ntyp,aref
c    &        ,conc,ncmp,ncmpx,iatm)
c     else
      call chklat(ibrav,a,coa,boa,alpha,beta,gamma,vc,r,g
     &        ,atmicv,atmicp,itype,natm,anclr,rmt,ntyp,aref
     &        ,conc,ncmp,ncmpx,iatm,asa,angl,fill)
c     endif
      write(*,1510) bravai(ibrav),a,coa,boa,alpha,beta,gamma
 1510 format(/'   lattice constant'
     &       /'   bravais=',a3,'   a=',f9.5,'   c/a=',f7.4,
     &        '   b/a=',f7.4/'   alpha=',f7.2,'   beta=',f7.2,
     &        '   gamma=',f7.2)
      write(*,'(a,f12.5,a)')'   unit cell volume=',vc*a**3,'(a.u.)'
      write(*,'(a,f6.1,a)')'   volume filling=',1d2*fill,'%'
      if(ifkey('tlt',brvtyp))
     &  write(*,'(a,3f11.5,a)')'   tilt angle='
     &                   ,angl(3),angl(2),angl(1),' (Euler angle)'
      write(*,'(/a)')'   primitive translation vectors (in unit of a)'
      write(*,'(a,3f9.5,a)')'   a=(',(r(i,1),i=1,3),')'
      write(*,'(a,3f9.5,a)')'   b=(',(r(i,2),i=1,3),')'
      write(*,'(a,3f9.5,a)')'   c=(',(r(i,3),i=1,3),')'
      write(*,'(/a)')'   reciprocal lattice vectors (in unit of 2*pi/a)'
      write(*,'(a,3f9.5,a)')'   ga=(',(g(i,1),i=1,3),')'
      write(*,'(a,3f9.5,a)')'   gb=(',(g(i,2),i=1,3),')'
      write(*,'(a,3f9.5,a)')'   gc=(',(g(i,3),i=1,3),')'
      write(*,'(/a)')'   type of site'
      ipair=0
      ji=0
      do 11 i=1,ntyp
      write(*,1210)type(i),rmt(i),field(i),lmxtyp(i)
 1210 format('   type=',a8,'  rmt=',f8.5,' field=',f7.3
     &      ,'   lmxtyp=',i3)
      jistrt=ji+1
      do 11 j=1,ncmp(i)
      ji=ji+1
c     --- find out pairs that are the same atomic species
c         in each type. In the LMD mode, those atoms are
c         treated to have magnetic moments of opossite
c         signs.
      do 310 jj=jistrt,ji-1
      if(abs(anclr(ji)-anclr(jj)) .le. 1d-6) then
      ipair=ipair+1
      if(ipair .gt. lmpmx) call errtrp(1,'spmain','ipair too large')
      lmpair(1,ipair)=jj
      lmpair(2,ipair)=ji
      endif
  310 continue
   11 write(*,1212)j,anclr(ji),conc(ji)
 1212 format(16x,'  component=',i2,'  anclr=',f5.0,'  conc=',f7.4)
      write(*,'(/a)') '   atoms in the unit cell'
      do 12 j=1,natm
      i=itype(j)
      rmtd(j)=rmt(i)
   12 write(*,1310)(atmicp(l,j),l=1,3),atmtyp(j)
 1310 format('   position=',3f13.8,'  type=',a8)
      if(ids .eq. 5) then
      call xtoken(token,*16)
      go to 17
   16 call errtrp(1,'spmain','fspin data not found')
   17 read(token,*)fspin
      write(*,'(1x)')
      write(*,'(a,f13.6)')'   fixed spin moment=',fspin
      endif
      write(*,'(1x)')
c
c
c     --- anclr data on file 05 are used --
c         (i.e. anclr on file 24 ignored)
      readin=.false.
c
c     --- 2nd record read in --
      if(record .eq. '2nd') then
      rewind(24)
c     --- branch if eof of file (unit=24) detected ---
      read(24,end=10,err=10)
      read(24,end=20,err=20)itr1,ncmp24,meshr,zdmy,corlvl,dr,xr,v1,ef
     &                     ,er0,ew,ez,edelt0,dmy
c     if(ncmp24 .ne. ncmpx) go to 20
      readin=.true.
      go to 30
   10 call errtrp(2,'spmain','eof detected; data generated')
      record='init'
      go to 30
   20 call errtrp(2,'spmain','eof detected; 1st data used')
      record='1st'
   30 continue
      endif
c
c     --- 1st record read in --
      if(record .eq. '1st') then
      rewind(24)
c     --- branch if eof of file (unit=24) detected ---
      read(24,end=40,err=40)itr1,ncmp24,meshr,zdmy,corlvl,dr,xr,v1,ef
     &                     ,er0,ew,ez,edelt0,dmy
c     if(ncmp24 .ne. ncmpx) go to 40
      readin=.true.
      go to 50
   40 call errtrp(2,'spmain','eof detected; data generated')
      record='init'
   50 continue
      endif
      if(record .eq. 'init') then
c     --- pre-assumed Fermi energy used if no data are available.
c         potential data will be generated later.
      if(asa) then
c     --- symmetry breaking field is not effective for ASA.
      ef(1)=ef0
      ef(2)=ef0
      else
      ef(1)=ef0+dex
      ef(2)=ef0-dex
      endif
      er0=0d0
      edelt0=0d0
      ew=0d0
      ez=0d0
      itr=0
      endif
      if(.not. readin .and. record .ne. 'init')
     &    call errtrp(1,'spmain','illegal token')
c
c     ef(1)=ef(1)-2.1d0
c     ef(2)=ef(2)-2.1d0
      efsp=max(ef(1),ef(2))+dble(e(mse,1))
      efif=min(ef(1),ef(2))+dble(e(mse,1))-ewidth
      er0=0d0
c     ew=0d0
c     ez=0d0
      edelt0=0d0
c------- modify ew and ez when they become unsafe
      emrlim=0.07d0
c     emrgn=0.2d0
      emrgn=0.7d0
      if(ew+ez .lt. efsp+emrlim .or.
     &   ew-ez .gt. efif-emrlim .or.
     &   ez .gt. ewidth+emrgn) then
      call errtrp(3,'spmain','new ew, ez generated')
      ew=5d-1*(efsp+efif)
      ez=5d-1*(efsp-efif)+emrgn
      endif
c------- modify ew and ez every time the calculation starts
c     e2u=efsp+2d-1
c     e1u=efif-2d-1
c     ew=5d-1*(e2u+e1u)
c     ez=5d-1*(e2u-e1u)
c------- fix ew and ez to some special values
c     write(*,'(1x,a)')'   Fixed ew and ez are used'
c     ew=-0.45563
c     ez=1.2
c---------------------------------------
      unit=2d0*pi/a
      refunt=2d0*pi/aref
      e1=(ew-ez)/unit**2
      e2=(ew+ez)/unit**2
      write(*,'(1x,a,f10.5,a,f10.5)')'  ew=',ew,'  ez=',ez
c
c     --- e3 is used to select the prefered set from the
c     --- reciprocal lattice vectors. 'refunt'
c     --- is 2d0*pi/aref, where 'aref' is a reference lattice
c     --- constant fixed  by use of tabulated atomic radii.
c
c     --- e4 is used to fix 'eta' parameter used in the Ewalt
c     --- sum for the structural Green function. Here, not only
c     --- energies in the positive domain but also in the negative
c     --- one should be take into account for a good convergence.
      e3=(ew+ez)/refunt**2
      e4=max(abs(ew+ez),abs(ew-ez))/refunt**2
      if(abs(a-aref)/aref .gt. 2d-1) then
      call errtrp(3,'spmain','aref voided')
      e3=(ew+ez)/unit**2
      e4=max(abs(ew+ez),abs(ew-ez))/unit**2
      endif
c     e4=e2
c
c
      nf=nfqlty(bzqlty,ibrav)
      call tchmta(ew,ez,elvl,tm,ng)
      call cgtabl(ncub,clks,lmax,last)
c     --- get optimal eta value and the length of the corresponding
c     --- largest lattice vector and the reciprocal lattice vector,
      call etaopt(e4,vc,eta,rmx,gmx)
c     --- then generate the lattice vectors used in the Ewalt sum
      nrpmx=nwk/3
      call genrpt(r,wk,nrpmx,nrpt,rmx,atmicp,natm)
      allocate(rpt(3*nrpt))
      call equarr(wk,rpt,3*nrpt)
c     --- obtain the vertex positions of the first BZ
c     --- and the reciprocal latttice vectors.
      call bzvrtx(g,gcornr,ncg)
      ngpmx=nwk/3
      call gengpt(g,wk,ngpmx,ngpt,gmx,gcornr,ncg)
      allocate(gpt(3*ngpt))
      call equarr(wk,gpt,3*ngpt)
c     --- generate the prefered and remainder set
      call gnpset(gcornr,ncg,e3+emxr,gpt,ngpt,np,nr)
c     --- Now calculate coefficients for the Madelung potential.
      call madlng(eta,vc,atmicp,natm,rpt(4),nrpt-1,gpt(4)
     &           ,ngpt-1,rmtd,amdlng,asa)
c
c     --- Generate rotation matrices of cubic or hexagonal symmetry
c         operations that are applied to vectors, i.e., real harmonics
c         of l=2.
      isys=1
      iopmx=0
      do 13 j=1,2
      if(j .eq. 1) then
      call cubmtr(uprot,2)
      else
      call hexmtr(uprot,2)
      endif
c     --- Construct rotation matrices applied to vectors, i.e.
c         the l=1 component of real harmonics (y,z,x).
      do 19 iop=1,24
      do 19 i=1,9
      protdm(i,iop,j)=uprot(i,2,iop)
   19 protdm(i,24+iop,j)=-uprot(i,2,iop)
c     --- Get symmetry operations compatible with a given
c         Bravais lattice.
      call brvsy2(j,isymdm(1,j),ls,r,g,protdm(1,1,j))
c     --- Which coordinate system is likely, cubic or hexagonal?
c         The one that has a larger number of symetry operations
c         is likely to be the suitable frame for the system.
      ioptot=0
      do 21 iop=1,24
   21 ioptot=ioptot+isymdm(iop,j)
      if(ioptot .gt. iopmx) then
      iopmx=ioptot
      isys=j
      endif
   13 continue
c     --- Now, whether the frame is cubic or hexagonal is fixed.
c         Generate rotation matrices of cubic or hexagonal symmetry
c         operations that are applied to all real harmonics.
      if(isys .eq. 1) then
c     --- likely a system with cubic frame
      call cubmtr(wk,mxl)
c     call cubmtr(urotat,mxl)
      else if(isys .eq. 2) then
c     --- likely a system with hexagonal frame
      call hexmtr(wk,mxl)
c     call hexmtr(urotat,mxl)
      else
      call errtrp(1,'spmain','illegal isys')
      endif
      if(mxl .eq. 1) then 
      do 25 i=1,24
   25 urotat(i,1,1)=(1d0,0d0)
      else 
      ii=0 
      do 24 iop=1,24
      do 24 l=1,mxl
      do 24 i=1,(2*mxl-1)**2
      ii=ii+1
   24 urotat(i,l,iop)=dcmplx(wk(ii),0d0)
c     --- copy information suitable for the system
      do 18 iop=1,48
      isymop(iop)=isymdm(iop,isys)
      do 18 i=1,9
   18 protat(i,iop)=protdm(i,iop,isys)
      endif
c
c     --- Now, we transform the rotation matrix u from the
c         real hermonics based to the spherical hermonics based
c         one if spin-orbit coupling is to be considered.
       call srtrns(uu,mxl)
      if(ls .eq. 1) call uutrns(urotat,uu,mxl,2)
c     --- Full Brillouin zone instead of irreducible one will be taken
c         when 'bzqlty' containes the key word 'fbz'.
      if(ifkey('fbz',bzqlty)) then 
c     --- if the symmetry operation is to be supressed the following
c         4 lines must be uncommented.
      do 33 i=1,24
      isymop(i)=0
   33 isymop(i+24)=0
      isymop(1)=1
      fbz=.true.
      else
      fbz=.false.
      endif
c     --- Give aribitrary numbers corresponding to the atomic
c         scattering factors.
      call getnum(1997,wk,ntyp)
      do 31 i=1,natm
   31 wk(ntyp+i)=wk(itype(i))
c     --- Get scattering geometric structure factor.
      allocate(gsf(ngpt))
      call stfact(atmicp,wk(ntyp+1),natm,gpt,gsf,ngpt)
c     --- Among the above, only those operations that are
c         also compatible with the atomic arrangement within
c         the unit-cell are retained.
      call atmrot(irotat,natm,g,atmicp,itype,protat,isymop
     &           ,wk(ntyp+natm+1),gsf,wk(ntyp+1)
     &           ,gpt,ngpt)
      deallocate(gsf)
c     call chktyp(irotat,itype,natm)
c     --- Based on the above information, the k-point mesh covering
c         the first BZ is constructed.
c     write(*,'(a,24i2/)')'   isymop=',isymop
      if(ids .eq. 4) then
c     --- determine read and write format for k-vectors and
c         Bloch spectrum function. The input may be such as ,e.g.,
c         'spc33', meaning that the input format type 3 and output
c         format 3 will be used.
      n=lftfil(go)
      if(n .eq. 3) then
c     --- default format is the (c) format.
      irdfmt=3
      iwrfmt=1
      else
c     --- 0<irdfmt<4 and 0<iwrfmt<4
      read(go(4:n),'(i2)')nbuff
      if(nbuff .le. 9) nbuff=nbuff*11
      irdfmt=nbuff/10
      iwrfmt=mod(nbuff,10)
      if(irdfmt .lt. 1 .or. irdfmt .gt. 3) irdfmt=3
      if(iwrfmt .lt. 1 .or. iwrfmt .gt. 3) iwrfmt=3
      endif
c     --- read-in k-point on which the spectrum function is
c         calculated.
c     irdfmt=3 corresponds to the c-format (see subroutine readk)
      nk3mx=nwk/3
      call readk(wk,coa,boa,nk3,nk3mx,kcrt,25,kblst,irdfmt)
      allocate(spctrl(mse*nk3))
      open(unit=27,file=spcfil(1),status='unknown',form='formatted')
      open(unit=28,file=spcfil(2),status='unknown',form='formatted')
      else
      nk3=0
      endif
      ibit=31
      if(4*nf**3/ibit+1 .gt. nwk-3*nk3)
     &     call errtrp(1,'spmain','nf too large')
      ist=3*nk3+1
      call bzmesh(nf,nk,g,r,protat,isymop,wk(ist),ls,fbz)
      nkx=nk+nk3
      allocate(tch((2*mxl-1)**2*ng*nkx*ndmx),pexf(np*nkx*ndmx)
     &        ,vkp(3*nkx),wtkp(nkx),prr(np*nkx),hh(np*(2*mxl-1)**2*nkx)
     &        ,korder(nkx),stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'spmain','allocation fails')
      call ftchvk(vkp,wtkp,g,nk,wk(ist),nf)
c     --- randomize the order of sequence of the k-point mesh.
      call rndkpt(korder,nk,iwk)
      nthread=numcor()
      if(nthread .gt. 16) nthread=nthread-1
      length=nk*mse/nthread+min(mod(nk*mse,nthread),1)
      nk1=nk
      nk=nk+nk3
      if(ids .eq. 4) then
      do 142 k=1,3*nk3
  142 vkp(3*nk1+k)=wk(k)
      do 143 k=nk1+1,nk1+nk3
  143 korder(k)=k
      endif
c
c     --- Calculate the coefficients used in constructing the
c         structural Green's function.
      call gtchst(eta,vc,nk,vkp,lmxtyp,np,ngpt,gpt,tch,prr,hh,e1,e2,ng
     &           ,(2*mxl-1)**2,rpt,atmicp,nrpt,pexf,iblk,natm,ndmx,inv
     &           ,nd,itype,lmxblk,itblk,its)
c     --- get the row or column positions where each atomic starts in
c         the packed form of matrices. 
      ip=1
      do 150 i=1,natm
      iatmp(i)=ip
  150 ip=ip+(lmxtyp(itype(i))+1)**2
c
      write(*,2100)last,np,np+nr,nrpt,nk,nd
 2100 format('   last=',i4,'   np=',i3,'   nt=',i4,'   nrpt=',i4,
     &       '   nk=',i5,'   nd=',i4/)
c
c     --- generate starting potential by atomic calculation
c     --- dex is the symmetry breaking field (void later
c     --- if magtyp='nmag')
      if(record .eq. 'init') then
      call gsdatp(wk,a,vc,atmicp,r,anclr,corlvl,ro,rmt,dr,xr
     &           ,iwk,itype,natm,ntyp,meshr
     &           ,ncmp,ncmpx,conc)
c     --- core levels are relative to eorg in the case of ASA.
c      if(asa) then
c      do 260 i=1,ntyp
c      do 260 j=1,18
c  260 corlvl(j,i,1)=corlvl(j,i,1)-eorg
c      endif
      call clrarr(ro(1,1,2),ncmpx*meshr)
      call equarr(corlvl(1,1,1),corlvl(1,1,2),18*ncmpx)
      call potenv(sdftyp,itype,ntyp,natm,anclr,ro,v1,q,exspl
     &           ,dr,xr,meshr,a,wk,amdlng,u,ncmp,ncmpx,conc)
c     write(*,'((1x,1p6e13.6))')(-5d-1*v1(k,1,1)*xr(k,1),k=1,60)
c     write(*,'((1x,1p6e13.6))')(v1(k,1,1),k=1,meshr,10)
      do 78 is=1,2
      do 78 i=1,ncmpx
   78 call finitn(v1(1,i,is),anclr(i),xr(1,i),meshr)
      readin=.true.
      endif
c
c     --- write on data as the 1st record --
      rewind(24)
      if(outtyp .eq. 'update') then
      write(*,1600)
 1600 format('   record 1 will be overlaied by input and'
     &      /'   record 2 will be replaced by new output.'/)
      write(24)itr1,ncmpx,meshr,zdmy,corlvl,dr,xr,v1
     &        ,ef,er0,ew,ez,edelt0,dmy
c     write(24)itr1,ncmpx,meshr,zdmy
c    &   ,((((corlvl(i,j,k),i=1,15),j=1,ncmpx)
c    &     ,k=1,2),dr,xr,v1,ef
      endif
c
c
c     --- fix the core state configuration
      call corcnf(anclr,config,corlvl,ncmpx,ef)
c
c     --- give up/down avaraged data for magtyp='nmag' cases
      if(magtyp .eq. 'nmag') then
      call avearr(v1(1,1,1),v1(1,1,2),ncmpx*meshr)
      call avearr(config(1,1,1),config(1,1,2),18*ncmpx)
      call avearr(corlvl(1,1,1),corlvl(1,1,2),18*ncmpx)
      call avearr(ef(1),ef(2),1)
c
c     --- swap up/down spin data if necessary
c     ---  (magtyp='rvrs' or '-mag') --
      else if(magtyp .eq. 'rvrs' .or. magtyp .eq. '-mag') then
      call swparr(v1(1,1,1),v1(1,1,2),ncmpx*meshr)
      call swparr(config(1,1,1),config(1,1,2),18*ncmpx)
      call swparr(corlvl(1,1,1),corlvl(1,1,2),18*ncmpx)
      call swparr(ef(1),ef(2),1)
c
c     --- give a kick which enforces the system magnetic
c     --- if magtyp='kick' is specified.
c     --- kick is not effective for asa.
      else if(magtyp .eq. 'kick' .and. (.not. asa)) then
      ef(1)=ef(1)+dex
      ef(2)=ef(2)-dex
      endif
c
c     --- get some information
      volnrm=0d0
      vint=vc/4d0/pi
      do 22 i=1,ntyp
   22 iwtyp(i)=0
      do 23 ia=1,natm
      i=itype(ia)
      iwtyp(i)=iwtyp(i)+1
      volin=rmt(i)**3
      volnrm=volnrm+volin
   23 vint=vint-volin/3d0
      vint=a**3*vint
c
c     --- check radial mesh if it is consistent with a. if not,
c         potential data are interpolated so as to fit new mesh.
      celvol=a**3*vc
      ji=0
      do 41 i=1,ntyp
      atvol=celvol*rmt(i)**3/volnrm
      rtin=a*rmt(i)
      ratm=(3d0*atvol/4d0/pi)**(1d0/3d0)
      do 41 j=1,ncmp(i)
      ji=ji+1
   41 call ckmesh(rtin,ratm,xr(1,ji),dr(1,ji),v1(1,ji,1)
     &            ,v1(1,ji,2),wk,meshr)
c
c     --- initializatione --
c     construct a title card.
      title='------'
      do 151 i=1,ntyp
      if(iwtyp(i) .ge. 1) then
      if(iwtyp(i) .eq. 1) then
      cwtyp=' '
      else
      write(cwtyp,'(a,i3)')'_',iwtyp(i)
      endif
      call chleng(title,len)
      title(len+1:)='@'//type(i)//cwtyp
c     title(len+1:)=type(i)//cwtyp
      idmy=lftfil(title)
      endif
  151 continue
      idmy=lftfil(title)
      title(1:6)=' '
      do 153 j=7,len+1
  153 if(title(j:j) .eq. '@') title(j:j)=' '
      call chleng(title,len)
      write(*,1500)title(1:len)
 1500 format('   ***** self-consistent iteration starts *****'/a)
      qlty=0d0
      pbeta=6d-1
      call setarr(tend,-30d0,10)
      do 51 is=1,2
      do 51 i=1,ncmpx
c     ---If finite nuclear size correction has not yet been made,
c        here might be a good chance to make it.
      if(-v1(1,i,is)*xr(1,i) .gt. 1.9d0*anclr(i)) then
      write(*,'(1x,a)')'   finite nuclear size correction made'
      call finitn(v1(1,i,is),anclr(i),xr(1,i),meshr)
      endif
   51 continue
c     --- check for empty lattice
c     call clrarr(v1,2*meshr*ncmpx)
c     call setarr(v1,-1d-2,2*meshr*ncmpx)
      call equarr(v1,v3,2*meshr*ncmpx)
c     write(*,'((1x,1p,6e13.6))')(v1(k,1,1),k=1,meshr-1,10)
c     write(*,'((1x,1p,6e13.6))')(v1(k,1,1),k=meshr-50,meshr)
c
      if(magtyp .eq. 'nmag') then
      nspin=1
      do 42 i=1,ntyp
   42 field(i)=0d0
      else
      nspin=2
      endif
c
c     ---decompose the field into uniform and staggerd components.
      ufield=0d0
      do 14 j=1,natm
      i=itype(j)
   14 ufield=ufield+field(i)
      ufield=ufield/dble(natm)
      do 15 i=1,ntyp
   15 field(i)=field(i)-ufield
c
c     --- v(meshr,ji,k) is used for the iteration of the offset ef.
c         similar replacement for v1 and v2 will take place later --
c         This, however, is excuted only when the fixed spin
c         moment procedure is not used.
      if(ids .ne. 5) then
      do 118 i=1,ncmpx
      do 118 k=1,2
      if(asa) then
      v3(meshr,i,k)=eorg
      else
      v3(meshr,i,k)=ef(k)
      endif
  118 continue
      endif
c
c     --- iteration cycle starts here ---
      call utimer(time,0)
      cnvq=1d3
      do 120 itr=1,maxitr
c     --- following 6 lines can be discarded for slow I/O systems.
      if(bckup .and. outtyp .eq. 'update') then
      write(24)itr1,ncmpx,meshr,anclr,corlvl,dr,xr,v1,ef,er,ew,ez,
     &         edelt,dmy
      rewind(24)
      read(24)
      endif
c     --- 
      if(asa .and. ids .ne. 5) then
      ef(1)=5d-1*(ef(1)+ef(2))
      ef(2)=ef(1)
      endif
      call clrari(match,2*18*ncmpx)
c     if(itr .eq. 10)write(*,'((1x,1p,6e13.6))')
c    &    (v1(k,1,1,1),k=1,meshr-1,10)
      itr0=itr
      itr1=itr1+1
      call clrarr(ro,2*ncmpx*meshr)
c
c     --- the uniform component of the field is taken into account
c         by shifting the up and down fermi levels,
      efs(1)=ef(1)+ufield
      efs(2)=ef(2)-ufield
      if(asa .and. itr .eq. 1 .and. abs(ufield) .lt. 1d-10 .and.
     &     (record .eq. 'init' .or. magtyp .eq. 'kick')) then
      efs(1)=ef(1)+dex
      efs(2)=ef(2)-dex
      endif
c     --- and the staggered components by shifting the potential zero.
      ji=0
      do 121 i=1,ntyp
      do 121 j=1,ncmp(i)
      ji=ji+1
      if(asa) then
      v1(meshr,ji,1)=eorg+field(i)
      v1(meshr,ji,2)=eorg-field(i)
      else
      v1(meshr,ji,1)=field(i)
      v1(meshr,ji,2)=-field(i)
      endif
  121 continue
c
c     --- check if the energy range is fully covered by ew and ez.
c         if not, stop the iteration and backup the data.
      ereset=.false.
      do 122 is=1,nspin
  122 if(efs(is)-eunder .lt. ew-ez
     &   .or. efs(is)-eunder+ewidth .gt. ew+ez) ereset=.true.
      if(ereset) then
      efsp=max(dble(e(mse,1)),dble(e(mse,1)))
      efif=min(dble(e(mse,1)),dble(e(mse,1)))-ewidth
      ew=5d-1*(efsp+efif)
      ez=5d-1*(efsp-efif)+emrgn
      call errtrp(3,'spmain','new ew and ez generated:')
      write(*,'(1x,a,f10.5,a,f10.5)')'  ew=',ew,'  ez=',ez
      e1=(ew-ez)/unit**2
      e2=(ew+ez)/unit**2
      call tchmta(ew,ez,elvl,tm,ng)
c     --- Recalculate the coefficients used in constructing the
c         structural Green's function.
      call gtchst(eta,vc,nk,vkp,lmxtyp,np,ngpt,gpt,tch,prr,hh,e1,e2,ng
     &           ,(2*mxl-1)**2,rpt,atmicp,nrpt,pexf,iblk,natm,ndmx,inv
     &           ,nd,itype,lmxblk,itblk,its)
      endif
c
c     --- perform kkr --
c     call clrarr(sm,2*ng*mmxl*ncmpx)
      do 70 is=1,nspin
      call drvmsh(efs(is),ewidth,edelt,ebtm(is),e(1,is),mse,ids)
      call cgnwt(e(1,is),wt(1,1,1,is),ng,mse,ew,ez)
c
c     --- start with single site properties.
      do 60 i=1,ncmpx
c     --- atomic t-matrices
      msg=itr .eq. maxitr
      if(ls .eq. 0) then
c     --- if ls=0, spin-orbit coupling is not considered.
      call phasea(v1(1,i,is),tchc(1,1,i),tchs(1,1,i)
     &         ,fcs(1,1,i),mxl,mxlcmp(i),elvl,tm,ew,ez,ng,dr(1,i)
     &         ,xr(1,i),meshr,rstr(1,i,is),wk,isr,msg)
      else
c     --- if ls is not 0, spin-orbit coupling is considered.
      call phaseb(v1(1,i,is),tchc(1,1,i),tchs(1,1,i)
     &         ,fcs(1,1,i),mxl,mxlcmp(i),elvl,tm,ew,ez,ng,dr(1,i)
     &         ,xr(1,i),meshr,rstr(1,i,is),wk,isr,is,msg)
      endif
c
c     --- polinomial fitting of 1/s, which mimics 1/s but has
c     --- no singularities.
      call fczero(cm(1,1,i),ew,ez,tchc(1,1,i),tchs(1,1,i)
     &            ,fcs(1,1,i),tm,elvl,ng,mmxl,mxlcmp(i)**2,xlim)
c
c     --- core states
      call cstate(v1(1,i,is),ro(1,i,is),rorg(1,i,is)
     &    ,corlvl(1,i,is),ancp(i,is),dr(1,i),xr(1,i)
     &    ,match(1,i,is),config(1,i,is),esic(i,is),meshr
     &    ,wk,isr,sdftyp,ebtm(is),asa,mxlcmp(i))
   60 continue
c
c     --- central part of kkr-cpa.
      call spckkr(wkc,ff,e(1,is),mse,cm,tchc,tchs,fcs
     &           ,ew,ez,mxl,ng,iblk,natm,ntyp,itype,vc
     &           ,isr,a,nk,np,tch,pexf,prr,hh,ncub,clks
     &           ,last,detl(1,is),wtkp,nd,tmt,phf,str,tc
     &           ,iwtyp,ids,length,gfree,ess,ncmp,ncmpx,conc
     &           ,tcpa(1,is),phase,korder,convrg,cnvq,urotat,uu
     &           ,irotat,isymop,ck,cj,iatm,ls,spctrl,nk3,lmxtyp
     &           ,mxlcmp,msiz,lmxblk,iatmp,itblk,its)
c     --- check cpa convergence
      cpacnv(is)=.true.
      do 52 k=1,mse
   52 cpacnv(is)=cpacnv(is) .and. convrg(k)
c
c     --- construct sm
      call gensm(e(1,is),wt(1,1,1,is),wkc,mmxl,mxlcmp,ncmpx,mse
     &          ,sm(1,1,1,is),ng)
c
c     --- partial and total density of states at the fermi level
      kk=mse
      do 65 i=1,ncmpx
      do 65 l=1,mxlcmp(i)**2
   65 dosef(l,i,is)=-dimag(wkc(l,i,kk))/pi
      total(is)=dimag(detl(kk,is))
      def(is)=dimag((detl(kk,is)-detl(kk-1,is))/(e(kk,is)-e(kk-1,is)))
c     write(*,*)def(is)
c
c-----------------------------------------------------------------------
c     --- print partial and the total DOS if required.
      if(ids .eq. 1 .or. ids .eq. 2 .or. ids .eq. 3) then
      estep=dble(e(2,is))-dble(e(1,is))
      do 69 i=1,ncmpx
      write(*,'(//1x,a,i2,a,i2)')'DOS of component',i
      do 69 k=1,kk
      xmd(i,k,1,is)=-dimag(wkc(2,i,k))/pi
cc    &    (-dimag(wkc(4,i,k))/pi)+(-dimag(wkc(2,i,k))/pi)
      xmd(i,k,2,is)=-dimag(wkc(4,i,k))/pi
cc    &    (-dimag(wkc(4,i,k))/pi)-(-dimag(wkc(2,i,k))/pi)
c  69 write(*,'(1x,f7.4,3x,9f8.4)') dble(e(k,is))-ef(is)
c    &      ,( -dimag(wkc(l,i,k))/pi,l=1,mxl**2)
      do 160 l=1,mxlcmp(i)
c     do 160 l=1,2
      do 160 m=1,2*(l-1)
  160 wkc(l**2,i,k)=wkc(l**2,i,k)+wkc(l**2-m,i,k)
c     wkc(5,i,k)=wkc(5,i,k)+wkc(6,i,k)+wkc(8,i,k)
c     wkc(7,i,k)=wkc(7,i,k)+wkc(9,i,k)
   69 write(*,'(1x,f7.4,3x,4f10.4)') dble(e(k,is))-ef(is)
     &      ,(-dimag(wkc(l**2,i,k))/pi,l=1,mxlcmp(i))
c    &      ,(-dimag(wkc(l,i,k))/pi,l=2,4)
c    &      , -dimag(wkc(5,i,k))/pi, -dimag(wkc(7,i,k))/pi
cc    &      , -dimag(wkc(2,i,k))/pi, -dimag(wkc(4,i,k))/pi
      write(*,'(//1x,a/(1x,f12.7,f13.5))')
     &      'total DOS',(dble(e(k,is))-estep/2d0-ef(is)
     &      ,dimag((detl(k,is)-detl(k-1,is))/(e(k,is)-e(k-1,is)))
     &        ,k=2,kk)
      write(*,'(//1x,a/(1x,f12.7,f13.5))')
     &      'integrated DOS',(dble(e(k,is))-ef(is)
     &       ,dimag(detl(k,is)),k=1,kk)
      else if(ids .eq. 4) then
c     iwrtfmt=3 corresponds to the c-format (see subroutine wrtspc)
      call wrtspc(spctrl,vkp,is,e,mse,kcrt,kblst,nk1,nk3,ef,unit,iwrfmt)
      endif
c-----------------------------------------------------------------------
c
   70 continue
c     --- paramgnetic case --
      if(magtyp .eq. 'nmag') then
      call equarr(rstr(1,1,1),rstr(1,1,2),meshr*mmxl*ng*ncmpx)
      call equarr(ro(1,1,1),ro(1,1,2),meshr*ncmpx)
      call equarr(rorg(1,1,1),rorg(1,1,2),20*ncmpx)
      call equarr(corlvl(1,1,1),corlvl(1,1,2),18*ncmpx)
      call equarr(config(1,1,1),config(1,1,2),18*ncmpx)
      call equari(match(1,1,1),match(1,1,2),18*ncmpx)
      call equarr(sm(1,1,1,1),sm(1,1,1,2),ng*mmxl*ncmpx)
      call equarc(detl(1,1),detl(1,2),mse)
      call equarr(dosef(1,1,1),dosef(1,1,2),mmxl*ncmpx)
      call equarc(e(1,1),e(1,2),mse)
      call equarc(wt(1,1,1,1),wt(1,1,1,2),ng*3*mse)
      call equarr(v1(1,1,1),v1(1,1,2),meshr*ncmpx)
      call equarr(ancp(1,1),ancp(1,2),ncmpx)
      call equarr(esic(1,1),esic(1,2),ncmpx)
      call equarc(tcpa(1,1),tcpa(1,2),mmxl**2*natm*mse)
      efs(2)=efs(1)
      ef(2)=ef(1)
      total(2)=total(1)
      def(2)=def(1)
      cpacnv(2)=cpacnv(1)
      endif
c
      if(ls .eq. 1 .and. ids .eq. 3) then
c     --- get matrix elements |<1s| r |p(e)>|**2 needed for the
c         calculation of the dipole transition.
      do 73 is=1,2
      do 73 i=1,ncmpx
c     --- K-edge XMD is calculated only when p valence states are available.
      if(mxlcmp(i) .ge. 2) then
c     write(*,'(1x,a,2i3)')'is,i=',is,i
      call dplini(v1(1,i,is),corlvl(1,i,is),rstr(1,i,is),tm,ng
     &           ,dr(1,i),xr(1,i),meshr,mmxl,wk,isr,match(1,i,is),asa)
      do 79 k=1,kk
      call dplmtx(dble(e(k,is)),ew,ez,dipole)
      do 79 l=1,2
   79 xmd(i,k,l,is)=1d5*dipole(l)*xmd(i,k,l,is)
      endif
   73 continue
      jbase=0
      do 75 i=1,ncmpx
c     do 75 i=1,ntyp
c     --- K-edge XMD is calculated only when p valence states are available.
c     if(lmxtyp(i) .ge. 2) then
      if(mxlcmp(i) .ge. 2) then
c     write(*,'(1x,a,i2)')'xmd of atom type',i
      write(*,'(1x,a,i2)')'xmd of component',i
      do 71 k=1,kk
   71 write(*,'(1x,f9.4,(2(2x,3f12.7)))')(dble(e(k,1))-ef(1))*13.6d0
     &   , xmd(i,k,2,1)+xmd(i,k,1,2)+xmd(i,k,1,1)+xmd(i,k,2,2)
     &    ,xmd(i,k,2,1)-xmd(i,k,1,1)
     &    ,xmd(i,k,1,2)-xmd(i,k,2,2)
c    &   ,(xmd(j,k,2,1)+xmd(j,k,1,2)+xmd(j,k,1,1)+xmd(j,k,2,2)
c    &    ,xmd(j,k,2,1)-xmd(j,k,1,1)
c    &    ,xmd(j,k,1,2)-xmd(j,k,2,2)
c    &    ,j=jbase+1,jbase+ncmp(i))
      endif
   75 jbase=jbase+ncmp(i)
c-----------------------------------------------------------------------
      endif
c
c     --- check charge neutrality, get the shift of the fermi level
c         and modify sm according to the shift.  ef, however, will
c         be shifted only after total energy calculation.
c     --- For fixed spin moment case, the fermi level is shifted
c         separately for each spin band.
      do 76 i=1,ncmpx
   76 anc(i)=ancp(i,1)+ancp(i,2)
      sold=total(1)-total(2)
      if(ids .eq. 5) then
c     --- ids=5 is the fixed spin moment case
      call fxspin(total,cnutr,anclr,anc,efs,dosef,def,sftef,sm
     &           ,ew,ez,itype,mmxl,mxlcmp,ng,natm,ntyp,ncmp,ncmpx
     &           ,conc,fspin)
      else
c     --- usual cases
      call neutrl(total,cnutr,anclr,anc,efs,dosef,def,sftef(1),sm
     &          ,ew,ez,itype,mmxl,mxlcmp,ng,natm,ntyp,ncmp,ncmpx,conc)
      sftef(2)=sftef(1)
      endif
c
c     --- construct new charge and spin densities --
c         variable f containes very important information though
c         it is not used in the main program presently. for
c         convenience of possible other type of calculations f is
c         stored in reserved area.
      do 72 is=1,2
      if(asa) then
      call chrasa(ro(1,1,is),rorg(1,1,is),tof(1,1,is)
     &           ,f(1,1,1,is),mmxl,mxlcmp,sm(1,1,1,is),tm,ng,ntyp,xr
     &           ,meshr,rstr(1,1,is),total(is),itype,natm
     &           ,ncmp,ncmpx,conc,lmxtyp)
      else
      call chrdnc(ro(1,1,is),rorg(1,1,is),tof(1,1,is)
     &           ,f(1,1,1,is),mmxl,mxlcmp,sm(1,1,1,is),tm,ng,ntyp,xr
     &           ,meshr,rstr(1,1,is),total(is),itype,natm
     &           ,ncmp,ncmpx,conc)
      endif
   72 continue
      call chkcnf(anclr,config,ncmpx,tof,mmxl,mxlcmp)
      call rrcnv(ro,ncmpx*meshr,ncmpx*meshr)
c
c     --- construct new potential --
      call potenv(sdftyp,itype,ntyp,natm,anclr,ro,v2,q,exspl
     &           ,dr,xr,meshr,a,wk,amdlng,u,ncmp,ncmpx,conc)
      do 77 is=1,2
      do 77 i=1,ncmpx
   77 call finitn(v2(1,i,is),anclr(i),xr(1,i),meshr)
c
c     --- calculate total energy --
      do 80 is=1,2
      kk=mse
      zgiven=detl(kk,is)+sftef(is)*(detl(kk,is)-detl(kk-1,is))
     &            /(e(kk,is)-e(kk-1,is))
      call banden(detl(1,is),e(1,is),wt(1,1,1,is),bnd2(is),kk,ng)
   80 bnd2(is)=bnd2(is)+dimag(e(kk,is)*zgiven)
      call totalw(ro,v1,wk,corlvl,bnd2,esic,meshr,xr,dr,natm,ntyp
     &           ,itype,u,te,config,ncmp,ncmpx,conc)
c
c     --- construct new input potential --
c         exchange splitting which is biased by the constant
c         shift of the muffin-tin potential is stored in the
c         v(meshr,i,j) and iterated together with the offset
c         potentials.
c     --- For the fixed moment case, the iteration for ef using
c         Tchebychev accelaration must be skipped.
      if(ids .ne. 5) then
      do 74 i=1,ncmpx
      if(asa) then
      v1(meshr,i,1)=eorg
      v1(meshr,i,2)=eorg
      v2(meshr,i,1)=eorg
      v2(meshr,i,2)=eorg
      else
      v1(meshr,i,1)=ef(1)
      v1(meshr,i,2)=ef(2)
      v2(meshr,i,1)=5d-1*(ef(1)+ef(2)+exspl)
      v2(meshr,i,2)=5d-1*(ef(1)+ef(2)-exspl)
      endif
   74 continue
      endif
c
c     --------- Tchebychev accelaration --------
      call erranc(v3,v1,v2,pmix,rate,pbeta,qlty,tend,rms,cnvq,dr,xr
     &           ,meshr,ncmpx)
      if(mod(itr,20) .eq. 0) then
      if(qlty .lt. -.95d0 .and. rate .lt. 1d0)
     &     pbeta=(1d0/rate-sqrt((1d0/rate)**2-1d0))**2
      if(qlty .gt.  0d0) pbeta=pbeta*.98d0
      endif
c
      ctotal=total(1)+total(2)
      stotal=total(1)-total(2)
      if(ids .eq. 5) stotal=sold
      cinter=ro(meshr,1,1)*vint
      sinter=ro(meshr,1,2)*vint
      status=' '
      if(cpacnv(1) .eqv. .false.)status(1:1)='*'
      if(cpacnv(2) .eqv. .false.)status(2:2)='*'
      write(*,1100)status,itr,cnutr,stotal,te,log10(cnvq)
 1100 format(1x,a2,'itr=',i3,' neu=',f10.4,'  moment=',f8.4
     &       ,'  te=',f16.8,'  err=',f7.3)
c
c     --- write on logging file if required (go='log') ---
      if(go .eq. 'log') then
      call uclock(cptime)
      open(26,file=logfil,form='formatted',
     &     status='unknown')
      do 99 irec=1,10000
   99 read(26,*,end=98)
      call errtrp(1,'spmain','too many records')
   98 backspace(26)
      write(26,1110)itr,cptime,cnutr,stotal,te,log10(cnvq)
 1110 format('it=',i3,' cp=',f7.0,' n=',f8.4,' m=',f8.4
     &       ,' e=',f16.8,' er=',f7.3)
      close(26)
      endif
c     --- procedure needed in the case of lmd
      nswap=0d0
      if(magtyp .eq. 'lmd') then
c     --- check magnetic moment of each pair
c         and swap data of the pair if they do not
c         have an oposite sign.
      exspl=5d-1*(ef(1)-ef(2))
      ef(1)=ef(1)-exspl
      ef(2)=ef(2)+exspl
      do 322 i=1,ncmpx
c     do 324 is=1,2
c     v1(meshr,i,is)=0d0
c     v2(meshr,i,is)=0d0
c 324 v3(meshr,i,is)=0d0
c     v1(meshr,i,is)=0d0
c 324 v2(meshr,i,is)=0d0
c 324 v3(meshr,i,is)=0d0
      do 326 k=1,18
      corlvl(k,i,1)=corlvl(k,i,1)-exspl
  326 corlvl(k,i,2)=corlvl(k,i,2)+exspl
c     do 322 k=1,meshr-1
      do 322 k=1,meshr
      v1(k,i,1)=v1(k,i,1)-exspl
      v2(k,i,1)=v2(k,i,1)-exspl
      v3(k,i,1)=v3(k,i,1)-exspl
      v1(k,i,2)=v1(k,i,2)+exspl
      v2(k,i,2)=v2(k,i,2)+exspl
  322 v3(k,i,2)=v3(k,i,2)+exspl
      do 320 i=1,ipair
      j1=lmpair(1,i)
      j2=lmpair(2,i)
      call jipinv(j1,ityp1,i1)
      call jipinv(j2,ityp2,i2)
      if(ityp1 .ne. ityp2) call errtrp(1,'spmain','lmd fails')
      spin1=0d0
      spin2=0d0
      do 330 j=1,(lmxtyp(ityp1)+1)**2
      spin1=spin1+tof(j,j1,1)-tof(j,j1,2)
  330 spin2=spin2+tof(j,j2,1)-tof(j,j2,2)
c     write(*,'(1x,a,2i2,a,2f12.5,a,f12.5)')'pair',j1,j2,
c    &  ' moment',spin1,spin2,' exspl',exspl
      if(spin1*spin2 .gt. 0d0) then
      nswap=nswap+1
      call swparr(v3(1,j2,1),v3(1,j2,2),meshr)
      call swparr(v2(1,j2,1),v2(1,j2,2),meshr)
      call swparr(v1(1,j2,1),v1(1,j2,2),meshr)
      call swparr(corlvl(1,j2,1),corlvl(1,j2,2),18)
      call swparr(config(1,j2,1),config(1,j2,2),18)
      endif
  320 continue
      endif
c
c     --- check convergence --
      if(cnvq .lt. tol) go to 130
c
c     --- modify ef --
      if(itr .lt. maxitr) then
      if(ids .eq. 5) then
c     --- iteration for ef is easy for the fixed spin moment case.
      redsft=5d-1
      do 112 is=1,2
  112 ef(is)=ef(is)+redsft*sftef(is)
      else
c     --- give a new fermi level. the averaging with respect to the
c         atom species is taken only formally.
      if(asa) then
      redsft=1d0
c     redsft=5d-1
      do 114 is=1,2
  114 ef(is)=ef(is)+redsft*sftef(is)
      else
      redsft=2d-1
      do 110 is=1,2
      ef(is)=redsft*sftef(is)
      do 110 i=1,ntyp
      do 110 j=1,ncmp(i)
      call jip(i,j,ji)
  110 ef(is)=ef(is)
     &   +conc(ji)*dble(iwtyp(i))*v1(meshr,ji,is)/dble(natm)
      endif
      endif
      endif
  120 continue
c
c     --- iteration loop ends here --
c
      write(*,1800)
 1800 format('   *** no convergence')
  130 call utimer(time,itr0)
      do 132 is=1,2
c      ef(is)=v3(meshr,i,is)
      do 132 i=1,ncmpx
      if(asa) then
      v3(meshr,i,is)=eorg
      else
      v3(meshr,i,is)=0d0
      endif
  132 continue
      write(*,2000)sdftyp,reltyp,pmix,title(1:len)
 2000 format('   sdftyp=',a5,'   reltyp=',a5,'   dmpc=',f8.5/a)
      write(*,1200)itr0,cnutr,ctotal,stotal,cinter,sinter
     &        ,((log10(rms(i,is)),i=1,ncmpx),is=1,2)
 1200 format('   itr=',i3,'  neu',f8.4,'  chr,spn',2f8.4
     &       ,'  intc,ints',2f8.4/'   rms err=',6f7.3
     &       /(11x,6f7.3))
      write(*,2500)ef,def,te
 2500 format('   ef=',2f12.7,'  def=',2f12.7
     &      /'   total energy=',f17.9)
      if(ids .eq. 5) write(*,'(1x,a,f12.7,a,f12.7)')
     &   '  fixed spin moment=',fspin,
     &   '  field used to fix the spin=',ef(1)-ef(2)-exspl
c
      write(*,'(///)')
c     stop
      do 140 i=1,ntyp
      do 140 j=1,ncmp(i)
      call jip(i,j,ji)
      call dsenum(type(i),anclr(ji),anc(ji),config(1,ji,1)
     &           ,config(1,ji,2),corlvl(1,ji,1),corlvl(1,ji,2)
     &           ,tof(1,ji,1),tof(1,ji,2),mxlcmp(ji),ls)
  140 call hypera(type(i),corlvl(1,ji,1),corlvl(1,ji,2),ew,ez
     &           ,v3(1,ji,1),v3(1,ji,2),dr(1,ji),xr(1,ji)
     &           ,rorg(1,ji,1),rorg(1,ji,2),meshr,config(1,ji,1)
     &           ,config(1,ji,2),isr,hhf(ji),anclr(ji))
c
      if(outtyp .eq. 'update') then
c
c     --- swap up/down spin data if stotal is negative
      if(stotal .lt. -1d-10) then
      call swparr(v3(1,1,1),v3(1,1,2),ncmpx*meshr)
      call swparr(config(1,1,1),config(1,1,2),18*ncmpx)
      call swparr(corlvl(1,1,1),corlvl(1,1,2),18*ncmpx)
      call swparr(ef(1),ef(2),1)
      endif
      write(24)itr1,ncmpx,meshr,anclr,corlvl,dr,xr,v3,ef,er,ew,ez,
     &         edelt,dmy
c     write(*,'(1x,a,i4)')'data written successfully  itr=',itr1
c     write(*,'((1x,1p,6e13.6))')(v3(k,1,1),k=1,meshr-1,10)
      endif
c    &   write(24)itr1,ncmpx,meshr,anclr
c    &   ,(((corlvl(i,j,l),i=1,15),j=1,ncmpx),l=1,2)
c    &   ,dr,xr,v3,ef
      rewind(24)
      write(*,'(///)')
      open(25,file=inffil,form='formatted',
     &     status='unknown')
c     --- in order to append new records, old ones are skipped
c         before writing on.
      do 95 irec=1,10000
   95 read(25,*,end=94)
      call errtrp(1,'spmain','too many records')
   94 backspace(25)
c--------  If necessary, special information is wrriten down on a
c          file in the following. This part should be castumized
c          for each case.
c
c-------- c/a dependence of the total energy is given
c     write(25,'(1x,2f12.2,f20.7)')a,coa,te
c
c-------- internal parameter dependence is given
c     write(25,'(1x,2f12.3,f20.7)')a,atmicp(3,3),te
c
c-------- relaxation dependence of te and Hhf are given
c     write(25,'(1x,2f12.3,f20.7,6f8.1)')atmicp(2,1),atmicp(1,2),te
c    &     ,(hhf(i),i=1,ncmpx)
c
c-------- a dependence of te and the moment are given
      write(25,'(1x,f12.4,f20.7,f10.5)')a,te,stotal
c
      close(25)
      deallocate(tch,pexf,vkp,wtkp,prr,hh,korder,gpt,rpt,stat=ierr)
      if(ierr .ne. 0) call errtrp(2,'spmain','deallocation fails')
      if(ids .eq. 4) then 
      deallocate(spctrl,stat=ierr)
      if(ierr .ne. 0) call errtrp(2,'spmain','deallocation fails')
      endif
      end
