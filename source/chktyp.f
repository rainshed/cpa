      subroutine chktyp(irotat,itype,natm)
c-----------------------------------------------------------------------
c     Given the table 'irotat', which describes the rotation of atoms
c     in the unit cell, and 'itype', which gives the type of each atom,
c     this routine tells if the information given by 'itype' is
c     compatible with 'irotat'. This is useful when the unit cell
c     contains many atom. It checks if the type assigned to each atom
c     is consistent with the crystal symmetry.
c     Coded by H. Akai, 1 Sept. 1999, Osaka.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer irotat(natm,24),itype(natm)
      logical*1,allocatable:: ini(:)
      allocate(ini(natm))
      do 10 i=1,natmx
   10 ini(i)=.true.
      do 20 i=1,natm
      if(ini(i)) then
      it=itype(i)
      do 30 j=1,natm
      if(itype(j) .eq. it) then
      do 40 l=1,24
      if(irotat(i,l) .eq. j) go to 30
   40 continue
      write(*,'(3x,a,i3,a,i3,a)')'atoms',i,' and',j,
     &  ' are assigned to the same type, but it obviously is wrong'
      call errtrp(1,'chktyp','check atomic position.')
      endif
   30 continue
      ini(j)=.false.
      endif
   20 continue
      deallocate(ini)
      end
