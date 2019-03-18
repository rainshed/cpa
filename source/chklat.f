      subroutine chklat(ibrav,a,coa,boa,alpha,beta,gamma,vc,r,g
     &         ,atmicx,atmicp,itype,natm,anclr,rmt,ntyp,aref
     &         ,conc,ncmp,ncmpx,iatm,asa,angl,fill)
c----------------------------------------------------------------------
c     This program check values of 'a' and 'rmt' and if they are
c     unrealistic (either very large or small) it trys to generate
c     some more reasonable values. For 'a', experimental or MJW data
c     are used depending on a>=1d+4 or a<=1d-4, respectively.
c     coded by H.Akai, April 1992, Osaka
c     revised by H.Akai, Feb. 1993, Osaka
c     revised by H.Akai, 12 Aug, 1995, Duisburg
c     revised by H.Akai, 14 Aug, 1997, Duisburg
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 anclr(ncmpx),conc(ncmpx),rmt(ntyp),r(3,3),g(3,3)
     &      ,atmicp(3,natm),sft(3),cartes(3),cn(3)
     &      ,xtilt(3,3),stilt(3,3),rtilt(3),angl(3)
      real*8,allocatable::dist(:,:),wsradi(:),atvol(:)
      character atmicx(3,natm)*24,buff*80
      integer itype(natm),iatm(ntyp),ncmp(ntyp)
      integer,allocatable::num(:)
      logical asa
      logical*1,allocatable::fixed(:)
      data zero/1d-10/, tiny/1d-3/, iop/2/
      allocate(dist(ntyp,ntyp),wsradi(ntyp),atvol(ntyp),num(ntyp)
     &        ,fixed(ntyp))
      pi=4d0*atan(1d0)
      rad=pi/180d0
c     rad=0d0
      if(iop .eq. 0) call errtrp(3,'chklat','iop=0 specified')
c     --- see how often each type appears and
c     --- asign atom number for each atom type.
      do 10 i=1,ntyp
   10 num(i)=0
      do 20 i=1,natm
      iatm(itype(i))=i
   20 num(itype(i))=num(itype(i))+1
c
c     --- choose type of bravais lattice if it is not given.
      if(ibrav .eq. 0) then
c     --- see which type dominates the structure.
      idom=0
      ji=0
      wmx=0
      do 30 i=1,ntyp
      do 30 j=1,ncmp(i)
      ji=ji+1
      w=conc(ji)*dble(num(i))
      if(w .gt. wmx) then
      idom=ji
      wmx=w
      endif
   30 continue
c     --- now fix the structure.
      ibrav=inqbrv(anclr(idom))
      endif
c
      call prmvec(ibrav,coa,boa,alpha,beta,gamma,vc,r,g,angl)
      call clrarr(xtilt,9)
      call clrarr(stilt,9)
      do 180 i=1,3
  180 stilt(i,i)=1d0
      xtilt(1,1)=1d0
      xtilt(2,2)=boa
      xtilt(3,3)=coa
      do 190 i=1,3
c     --- tilt the cartesian coordinate
      call vrotat(rtilt,stilt(1,i),rad*angl(3),rad*angl(2),rad*angl(1))
      do 200 j=1,3
  200 stilt(j,i)=rtilt(j)
c     --- tilt the cartesian coordinate scaled by c/a and b/a
      call vrotat(rtilt,xtilt(1,i),rad*angl(3),rad*angl(2),rad*angl(1))
      do 190 j=1,3
  190 xtilt(j,i)=rtilt(j)
c     --- a recommended value of atomic volume.
c     idata=2 for mjw values, 1 for experimental ones.
c     get average atomic volumes, AS radii, and the unit cell volume.
      idata=2
      if(a .gt. 0.99d6) idata=1
      uvol=0d0
      ji=0
      do 40 i=1,ntyp
      atvol(i)=0d0
      do 50 j=1,ncmp(i)
      ji=ji+1
   50 atvol(i)=atvol(i)+qvolum(anclr(ji),idata)*conc(ji)
      uvol=uvol+dble(num(i))*atvol(i)
      if(atvol(i) .gt. zero) then
      wsradi(i)=(atvol(i)*3d0/4d0/pi)**(1d0/3d0)
      else
      wsradi(i)=zero
      endif
   40 continue
c     --- give lattice constant if it is not given.
      aref=(uvol/vc)**(1d0/3d0)
      if(a .lt. 1.01d-6 .or. a .gt. 0.99d6) a=aref
c
c     --- convert the data refered to a crystal axis cordinate
c     --- to those refered to the cartesian cordinate. data forms
c     --- such as  1d0, 0.5d0a, b, 3c, x, 0.333333y are allowed.
      ia=ichar('a')
      iz=ichar('z')
      do 52 j=1,natm
      do 53 i=1,3
   53 cartes(i)=0d0
      do 54 i=1,3
      buff=atmicx(i,j)
c     ---get the position of the last non-blank character
      call chleng(buff,l)
      ibuff=ichar(buff(l:l))
c     ---check if the last character is a, b, c, x, y, or z.
      if(ibuff .ge. ia .and. ibuff .le. iz) then
c     ---the data is a multiple of the primitive vectors a, b, c
c     ---or the vectors x=(1,0,0), (b/a)y=(0,b/a,0), and
c        (c/a)z=(0,0,c/a).
      buff(l:l)=' '
      if(l .eq. 1) then
c     ---containing no numeric data. unity is assumed.
      atmicp(i,j)=1d0
      else
      atmicp(i,j)=redata(buff)
      endif
c
      if(ibuff .le. ia+2) then
c     ---it is a, b, or c.
      ibuff=ibuff-ia+1
      do 56 l=1,3
   56 cartes(l)=cartes(l)+atmicp(i,j)*r(l,ibuff)
      else if(ibuff .ge. iz-2) then
c     ---it is x, y, or z.
      ibuff=ibuff-iz+3
      do 57 l=1,3
   57 cartes(l)=cartes(l)+atmicp(i,j)*xtilt(l,ibuff)
      else
c     ---it is neither of them.
      call errtrp(0,'chklat','illegal character appears')
      write(*,'(1x,a)') buff(l:l)
      stop
      endif
c
      else
c     ---the data is a number, meaning multiple of (1,0,0), etc.
      atmicp(i,j)=redata(buff)
      do 58 l=1,3
   58 cartes(l)=cartes(l)+atmicp(i,j)*stilt(l,i)
      endif
   54 continue
      do 52 i=1,3
   52 atmicp(i,j)=cartes(i)
c     --- Relocate all atoms inside the primitive unit cell.
      do 220 j=1,natm
      do 210 i=1,3     
      cn(i)=0d0
      do 210 n=1,3
  210 cn(i)=cn(i)+atmicp(n,j)*g(n,i)
      do 220 i=1,3
      icn=int(cn(i)+1d3+1d-3)-1000
      do 220 n=1,3
  220 atmicp(n,j)=atmicp(n,j)-icn*r(n,i)
c
c     --- The muffin-tin radii are fixed in the following blocks
      do 60 i=1,ntyp
      if(rmt(i) .gt. zero) then
      fixed(i)=.true.
      else
      fixed(i)=.false.
      endif
   60 continue
      do 70 it=1,ntyp
      i=iatm(it)
c     --- sarch the shortest distance for each type of the atom.
      do 80 jt=1,ntyp
   80 dist(it,jt)=1d10
      do 70 ix=-1,1
      do 70 iy=-1,1
      do 70 iz=-1,1
c     --- give translation to generate atoms of different cells
      do 90 l=1,3
   90 sft(l)=dble(ix)*r(l,1)+dble(iy)*r(l,2)+dble(iz)*r(l,3)
      do 70 j=1,natm
      jt=itype(j)
c     --- of course, we should exclude the target atom itself.
      if((.not. (ix .eq. 0 .and. iy .eq. 0 .and. iz .eq. 0
     &   .and. j .eq. i)) .and. jt .ge. it) then
      dd=0d0
      do 100 l=1,3
  100 dd=dd+(atmicp(l,j)+sft(l)-atmicp(l,i))**2
      dd=sqrt(dd)
c     write(*,'(1x,a,3i3)')'type,atom,type target atom=',it,j,jt
c     write(*,'(1x,a,3i3,3f10.5)')'ix,iy,iz,sft=',ix,iy,iz,sft
c     write(*,'(1x,a,1p2e13.5)')'distance,dist=',dd,dist(it,jt)
      if(dd .lt. dist(it,jt)) then
      dist(it,jt)=dd
c     --- dist gives the minimum distance between it and jt-th
c     --- type of atoms
      endif
      endif
   70 continue
c     write(*,'(1x,2i3,1p,e13.6)')
c    & ((it,jt,dist(it,jt),it=1,ntyp),jt=1,ntyp)
c     --- check consistency of the muffin-tin radii, giving warnings
c     --- and redefining them if they are not consistent.
c     --- If iop=0 I do not care about such conflictions.
      if(iop .ne. 0) then
      do 130 it=1,ntyp
      do 130 jt=it,ntyp
      if(rmt(it)+rmt(jt) .gt. dist(it,jt)+zero) then
        write(*,'(3x,a,i2,a,f10.5,a,i2,a,f10.5)')
     &   'rmt(',it,')=',rmt(it),'  rmt(',jt,')=',rmt(jt)
        call errtrp(2,'chklat','given rmt''s conflict; reduced')
        red=dist(it,jt)/(rmt(it)+rmt(jt))
        if(iop .eq. 1) then
c     --- for iop=1, rmt's which conflict will be reduced with
c     --- keeping other rmt's unchanged.
c     --- take care when it and jt are equivalent.
          rmtit=red*rmt(it)
          rmt(jt)=red*rmt(jt)
          rmt(it)=rmtit
        else if(iop .eq. 2) then
c     --- for iop=2, all rmt's are reduced with keeping the
c     --- ratios between them constant.
          do 160 i=1,ntyp
  160     rmt(i)=red*rmt(i)
        else
          call errtrp(1,'chklat','illegal iop specified')
        endif
      endif
  130 continue
      endif
c     --- then we try to determine rmt's starting from the tightest
c     --- case and then proceed step by step.
      do 110 itry=1,ntyp+1
      redm=1d10
      do 120 it=1,ntyp
      do 120 jt=it,ntyp
      if(fixed(it)) then
        if(fixed(jt)) then
c     --- both rmt(it) and rmt(jt) are already fixed.
          go to 120
        else
c     --- rmt(it) is already fixed but rmt(jt) is not.
          red=(dist(it,jt)-rmt(it))/wsradi(jt)
      if(red .lt. 1d-3) red=1d10
        endif
      else
        if(fixed(jt)) then
c     --- rmt(jt) is already fixed but rmt(it) is not.
          red=(dist(it,jt)-rmt(jt))/wsradi(it)
      if(red .lt. 1d-3) red=1d10
        else
c     --- neither of rmt(it) and rmt(jt) are fixed.
          red=dist(it,jt)/(wsradi(it)+wsradi(jt))
        endif
      endif
      if(red .lt. redm) then
      redm=red
      itm=it
      jtm=jt
      endif
  120 continue
c     if(redm .gt. 9d9) return
      if(redm .gt. 9d9) go to 150
      if(redm .lt. tiny)
     &   call errtrp(1,'chklat','inadequate rmt specified')
      if(.not. fixed(itm)) then
      rmt(itm)=wsradi(itm)*redm
      fixed(itm)=.true.
      endif
      if(.not. fixed(jtm)) then
      rmt(jtm)=wsradi(jtm)*redm
      fixed(jtm)=.true.
      endif
c     write(*,*)'itm,jtm',itm,jtm,'rmt',rmt(itm),rmt(jtm)
  110 continue
      call errtrp(3,'chklat','fails in determining rmt')
  150 vtin=0d0
      do 140 i=1,ntyp
  140 vtin=vtin+dble(num(i))*rmt(i)**3
      fill=(vtin*4d0*pi/3d0)/vc
      if(asa) then
      if(iop .ne. 0) then
      red=fill**(-1d0/3d0)
      fill=1d0
c     write(*,*)'vc=',vc
c     write(*,*)'sum of vtin=',vtin*4d0*pi/3d0
      do 170 i=1,ntyp
  170 rmt(i)=rmt(i)*red
      endif
c 170 rmt(i)=wsradi(i)/a
c 170 write(*,'(3x,a,i2,a,f12.7)')'rmt(',i,')=',rmt(i)
      endif
c     write(*,'(3x,a,f12.7)')'volume filling=',fill
      deallocate(dist,wsradi,atvol,num,fixed)
      end
