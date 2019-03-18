      subroutine readin(go,file,brvtyp,a,coa,boa,alpha,beta,gamma
     &             ,edelt,ewidth,reltyp,sdftyp,magtyp,record,outtyp
     &             ,bzqlty,maxitr,pmxtyp,ntyp,type,ncmp,rmt,field,lmxtyp
     &             ,anclr,conc,natm,atmicx,atmtyp,r,angl,*)
c-----------------------------------------------------------------------
c     --- KKR-CPA version ---
c
c     The following input data are read in.
c
c     go:      process is to be executed or not (go/ngo/dos/dsp).
c     file:    file name used.
c     brvtyp:  type of the bravais lattice (fcc/bcc/hcp/sc/bct/st
c              /fco/bco/bso/so/bsm/sm/trc/rhb/fct/trg/prv/aux).
c              if brvtyp is "prv" or "aux", it means that the primitive
c              vectors will be given as the subsequent data.
c              In this case, the input data for c/a, b/a,
c              alpha, beta, and gamma should not appear but
c              just be skipped.
c     r:       Three primitive vectors. These vectors should
c              be input only when brvtyp="prv" or "aux" is specified.
c     a:       lattice constants in atomic unit along a axis.
c     c/a:     c/a ratio.
c     b/a:     b/a ratio.
c     alpha:   angle between b and c axis in degree.
c     beta:    angle between c and a axis in degree.
c     gamma:   angle between a and b axis in degree.
c     edelt:   small imaginary part attached to the Fermi energy.
c     ewidth:   width of the energy window covered by the energy contour.
c     reltyp:  type of relativistic treatment (nrl/sra/nrlls/srals).
c     sdftyp:  type of the parametrization of xc energy (vbh/mjw/vwn
c              /lmmjw/lmvbh/pymjw/pyvbh/pyvwn). Switch 'asa' can be
c              added, e.g., 'mjwasa' means that the atomic sphere
c              approximation (ASA) is to be exploited.
c     magtyp:  magnetic state (mag/nmag/rvrs/-mag/kick) rvrs and -mag
c              are equvalent. kick is used to kick off the system
c              to be magnetic.
c     record:  record used as input potential (init/1st/2nd).
c     outtyp:  whether to update potential file or not (update/quit).
c     bzqlty:  quality of BZ mesh (t/l/m/h/u or any integer number).
c     maxitr:  maximum number of the iteration.
c     pmxtyp:  combination of mixing parameter and mixing type.
c              for example 0.03Tchebyshev, 0.03ch, 0.03Tch, etc.
c              also 0.03Broyden, 0.03br, 0.03bry, etc. are accepted.
c              The mixing type can be omitted, like '0.03'. In this
c              case a defualt mixing type (defined by the callin
c              routine) is assumed.
c
c     ntyp:    number of inequivalent sites called 'type'.
c     type:    name, starting by alphabet, specifing the type.
c     ncmp:    number of component atoms on this site
c     rmt:     muffin-tin radius for this site.
c     field:   external local magnetic field on this site.
c     lmxtyp:  maximum l value considered for this site.
c     anclr:   nuclear charge of the component atom.
c     conc:    relative probability of the occupation of the component
c              atom at the site.
c                     .........
c              above two data are repeated untill all the component
c              atoms on this site are given.
c              Then the set of data (type, ncmp, rmt, field,
c              anclr, conc, anclr, conc,...) must be repeated
c              for all the types.
c
c     natm:    number of atoms in the unit cell.
c     atmicx:  atomic position (x,y,z) in the unit cell (in unit of a).
c              or (x*v1, y*v2, z*v3), where v1,v2,or v3 are a,b, or c.
c              here a,b, and c are the primitive vectors. e.g.
c              0.5 0.5 0.5  or 0.5x 0.5y 0.5z means (0.5,0.5,0.5), but
c              0.5a 0.5b 0.5c means 0.5a+0.5b+0.5c. also
c              0.5 0.5 0.5c and 0.5b 0.5a 0.5c, etc. are allowed.
c     atmtyp:  name of the type given above.
c                     .........
c              Above two data must be repeated until all the atoms
c              in the unit cell are completed.
c
c     Return 1 is excuted at eof.
c     coded by H.Akai, April 1992, Osaka
c     KKR-CPA version coded by H.Akai, 7 Sep. 1996. Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 anclr(*),rmt(*),field(*),conc(*),r(3,3),angl(3)
      integer ncmp(*),lmxtyp(*)
      character token*256,go*(*),file*(*),brvtyp*(*),reltyp*(*)
     &         ,sdftyp*(*),magtyp*(*),outtyp*(*),type(*)*(*)
     &         ,atmtyp(*)*(*),bzqlty*(*),record*(*),lwcase*256
     &         ,atmicx(3,*)*(*),pmxtyp*(*)
      logical ifkey
c
c     first i give the default values for the input data
      go         =   'go'
      file       =   'tmpfile'
      brvtyp     =   ' '
      a          =   0d0
      coa        =   0d0
      boa        =   0d0
      alpha      =   0d0
      beta       =   0d0
      gamma      =   0d0
      edelt      =   1d-3
      ewidth     =   1.2d0
      reltyp     =   'nrl'
      sdftyp     =   'mjw'
      magtyp     =   'mag'
      record     =   '2nd'
      outtyp     =   'update'
      bzqlty     =   't'
      maxitr     =   40
      pmxtyp     =   '0.024d0ch'
      ntyp       =   1
      ncmp(1)    =   1
      type(1)    =   'undef'
      rmt(1)     =   0d0
      field(1)   =   0d0
      lmxtyp(1)  =   2
      anclr(1)   =   0d0
      conc(1)    =   1d0
      natm       =   1
      atmicx(1,1)=   '0d0'
      atmicx(2,1)=   '0d0'
      atmicx(3,1)=   '0d0'
      atmtyp(1)  =   'undef'
      r(1,1)     =    1d0
      r(2,1)     =    0d0
      r(3,1)     =    0d0
      r(1,2)     =    0d0
      r(2,2)     =    1d0
      r(3,2)     =    0d0
      r(1,3)     =    0d0
      r(2,3)     =    0d0
      r(3,3)     =    1d0
      angl(1)    =    0d0
      angl(2)    =    0d0
      angl(3)    =    0d0
c
c     then i try to get all the values from console.
      call xtoken(token,*10)
      if(token .ne. ' ') go=lwcase(token)
      if(go .eq. 'end' .or. go .eq. 'eof') return 1
      call xtoken(token,*10)
      if(token .ne. ' ') file=token
      call xtoken(token,*10)
      if(token .ne. ' ') brvtyp=lwcase(token)
      if(ifkey('prv',brvtyp) .or. ifkey('aux',brvtyp)) then
      do 20 j=1,3
      do 20 i=1,3
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*)r(i,j)
   20 continue
      endif
      if(ifkey('tlt',brvtyp)) then
      do 30 i=1,3
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*)angl(4-i)
   30 continue
      endif
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) a
      if(.not. (ifkey('prv',brvtyp) .or. ifkey('aux', brvtyp))) then
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) coa
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) boa
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) alpha
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) beta
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) gamma
      endif
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) edelt
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) ewidth
      call xtoken(token,*10)
      if(token .ne. ' ') reltyp=lwcase(token)
      call xtoken(token,*10)
      if(token .ne. ' ') sdftyp=lwcase(token)
      call xtoken(token,*10)
      if(token .ne. ' ') magtyp=lwcase(token)
      call xtoken(token,*10)
      if(token .ne. ' ') record=lwcase(token)
      call xtoken(token,*10)
      if(token .ne. ' ') outtyp=lwcase(token)
      call xtoken(token,*10)
      if(token .ne. ' ') bzqlty=lwcase(token)
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) maxitr
      call xtoken(token,*10)
      if(token .ne. ' ') pmxtyp=lwcase(token)
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) ntyp
      icount=0
c
c     loop over distinct types of site
      do 40 i=1,ntyp
      call xtoken(token,*10)
      type(i)='undef'
      if(token .ne. ' ') type(i)=token
      call xtoken(token,*10)
      ncmp(i)=1
      if(token .ne. ' ') read(token,*) ncmp(i)
      call xtoken(token,*10)
      rmt(i)=0d0
      if(token .ne. ' ') read(token,*) rmt(i)
      call xtoken(token,*10)
      field(i)=0d0
      if(token .ne. ' ') read(token,*) field(i)
      call xtoken(token,*10)
      lmxtyp(i)=2
      if(token .ne. ' ') read(token,*) lmxtyp(i)
c
c     each site is occupied by ncmp(i) component atoms
      do 40 j=1,ncmp(i)
      icount=icount+1
      call xtoken(token,*10)
      anclr(icount)=0d0
      if(token .ne. ' ') read(token,*) anclr(icount)
      call xtoken(token,*10)
      conc(icount)=1d0
   40 if(token .ne. ' ') read(token,*) conc(icount)
c
      call xtoken(token,*10)
      if(token .ne. ' ') read(token,*) natm
c
c     loop over atoms in the unit cell
      do 50 i=1,natm
      do 60 j=1,3
      call xtoken(token,*10)
      atmicx(j,i)='0d0'
   60 if(token .ne. ' ') atmicx(j,i)=lwcase(token)
      call xtoken(token,*10)
      atmtyp(i)='undef'
   50 if(token .ne. ' ') atmtyp(i)=token
      return
   10 return 1
      end
