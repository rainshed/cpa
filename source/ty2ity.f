      subroutine ty2ity(atmtyp,natm,type,ntyp,itype)
c----------------------------------------------------------------------
c     Convert 'atmtyp' data to 'itype' data which is integer.
c     This establishes a mapping betwee the atomic sites and the atom
c     types. For example, suppose i denotes an atomic site in the
c     unit-cell, then the type of the atom occupying that site is
c     given by j=itype(i).
c     coded by H.Akai
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character atmtyp(natm)*(*),type(ntyp)*(*)
      integer itype(natm)
      do 10 i=1,natm
      do 20 j=1,ntyp
      if(atmtyp(i) .eq. type(j)) then
      itype(i)=j
      go to 10
      endif
   20 continue
      write(6,'(1x,a,a12)')'  atom type =',atmtyp(i)
      call errtrp(1,'ty2ity','the above type not defined')
   10 continue
      return
      end
