      subroutine atmrot(irotat,natm,g,atmicp,itype,protat,isymop,wk
     &                 ,gsf,ftype,gpt,ngpt)
c-----------------------------------------------------------------------
c     This program construct a table describing the mapping, which
c     occurs when the symmetry operation designated by iop is applied,
c     among atoms in the unit cell. The measnin of the table is, e.g.,
c     if lineary aligned five atoms 1,2,3,4,5 change their position
c     such that it looks like 4,5,1,2,3 after the operation, the vector
c     'irotat' carries values (4,5,1,2,3). If the operation iop is
c     allowed, isymop(iop) is set to 1, otherwise 0.
c
c     Note: Since the procedure performed here is a bit complicated, the
c     style of the program seems to be very ellegant, namely a lot of
c     "go to" statement are used. This however is the easiest way to
c     realize the present algorithm, which might be changed in future.
c
c     Coded by H.Akai, 23 Dec. 96, Osaka
c     Revised on 7 Aug. 1997, Duisburg
c     Revised on 8 May  2003, Osaka
c     Revised on 1 May 2011 and 23 July 2014.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 protat(9,48),atmicp(3,natm),g(3,3),wk(3,natm),c(3)
     &      ,gpt(3,ngpt),ftype(natm),dsplmt(3),gsf(ngpt)
      integer itype(natm),irotat(natm,48),isymop(48)
      logical eqvlat,str
      data zero/1d-8/,iprint/1/
c
c     write(*,'(/1x,a)')'iop g/u    irotat'
c     do 200 iop=2,48
c 200 isymop(iop)=0
      call clrari(irotat,natm*48)
      do 30 i=1,natm
   30 irotat(i,1)=i
      isymop(1)=1
      do 10 iop=2,48
c     --- Consider only the symetry operations compatible with 
c         those of the Bravais lattice. 
      if(isymop(iop) .ne. 0) then
c     --- First, reset "isymop" that will be later set when the local
c         symmetry is compatible with the symmetry operation "iop".
      isymop(iop)=0
c     --- Rotate the cartesian cordinates of all the atoms.
      call rotatm(atmicp,wk,protat(1,iop),natm)
c     --- Check if the rotated crystal produces the same
c         atomic structure factor.
      str=eqvlat(wk,natm,gpt,ftype,gsf,ngpt)
c     write(*,*)iop,str
c     --- If it is not the case, skip the remaining procedure.
      if(str) then
c     write(*,*)'iop,str=',iop,str
c     write(*,'((1x,9f7.3))')((wk(i,j),i=1,3),j=1,natm)
c     --- First, fix the shift of the coordinate for a trial pair.
c         The first atom of the this pair can always be rotated atom 1.
      do 60 id=1,natm
      call clrari(irotat(1,iop),natm)
c     --- The pair members should be of the same type.
      if(itype(id) .eq. itype(1)) then
c     --- try (id,1) pair.
      irotat(id,iop)=1
      do 70 i=1,3
c     --- A displacement is assumed.
   70 dsplmt(i)=wk(i,1)-atmicp(i,id)
c     --- Now a displacement compatible to (id,1) pair has been set.   
c     write(*,'(1x,a,i2,a,3f10.5)')'target=',id,'   dsp=',dsplmt
      do 80 ib=2,natm
c     --- For the rotated ib-th atom, check if any unrotated atom comes
c         to the same position.
      do 90 ia=1,natm
c     --- (ia,ib) pair is examined only when ia has not yet been paired
c          and also its type is the same as that of ib.
      if(irotat(ia,iop) .eq. 0  .and. itype(ia) .eq. itype(ib)) then
c     --- Check the shift for the (ia, ib) pair.
      do 100 j=1,3
  100 c(j)=0d0
      do 110 i=1,3
      d=wk(i,ib)-atmicp(i,ia)-dsplmt(i)
c     --- inner product d*g gives coeffecients c's that satisfy
c         d=c(1)*r(1)+c(2)*r(2)+c(3)*r(3), where r's are the primitive
c         lattice vectors.
      do 110 j=1,3
  110 c(j)=c(j)+d*g(i,j)
c     write(*,'(1x,a,2i3,a,3f10.5)')'pair=',ia,ib,'  shift=',c
      do 120 i=1,3
c     --- check if all c's are integer. if not, ia and ib-th atom do not
c         sit at the same place.
      if(abs(c(i)-dnint(c(i))) .gt. 1d-3) go to 90
  120 continue
      irotat(ia,iop)=ib
      go to 80
      endif
   90 continue
      go to 60
   80 continue
      go to 130
      endif
   60 continue
      go to 20
c     --- since all the atoms are found to have their partner,
c         "isymop" is now set. This is done either direct or indirect
c         rotation compatible with the symmetry.
  130 isymop(iop)=1
c     write(*,'(1x,2i3,3x,5(5i3,1x))')
c    &      iop,isymop(iop),(irotat(i,iop),i=1,natm)
      go to 10
   20 continue
c     call errtrp(3,'atmrot','eqvlat detected no structure differences')
      isymop(iop)=0
      endif
      endif
   10 continue
c
c    --- optional output
      if(iprint .eq. 1) then 
c    --- For cubic operation, protat(3,13) is not zero
      if(abs(protat(3,13)) .lt. zero) then
      write(*,'(/t4,a,t12,a,t14,a,t17,a,t20,a,t23,a,t27,a,t35,a
c    &      /t4,40(''-'')
     &      )')
     &       'symop','E','C6','C3','C2','C3-','C6-','C2''*6'
      write(*,'(2(t6,a,t11,i2,t14,i2,t17,i2,t20,i2,t23,i2,t27,
     &    i2,t31,6i2/))')'g',(isymop(i),i=1,12),'u',(isymop(i),i=25,36)
      else    
      write(*,'(/t4,a,t12,a,t15,a,t22,a,t28,a,t37,a,t45,a,t57,a
c    &      /t4,62(''-'')
     &     )')
     &       'symop','E','C4*3','C2*3','C4^3*3','C3*4','C3^2*4','C2''*6'
      write(*,'(2(t6,a,t11,i2,t14,3i2,t21,3i2,t28,3i2,t35,4i2,t44,
     &   4i2,t53,6i2/))')'g',(isymop(i),i=1,24),'u',(isymop(i),i=25,48)
      endif   
      endif   
      end
