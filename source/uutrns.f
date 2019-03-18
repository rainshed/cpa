      subroutine uutrns(urotat,u,mxl,il)
c---------------------------------------------------------------------
c     Transformation of rotation matrix from spherical harmonics
c     representation to real harmonics representation
c     il=1: spherical to real r=u*s*u^+
c     il=2: real to spherical s=u^+*r*u
c     Coded by H. Akai, 9 April 2014, Tokyo
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 urotat((2*mxl-1)**2,mxl,24),u(2*mxl-1,2*mxl-1,2)
      complex*16,allocatable::wk(:,:)
      data zero/1d-20/
      allocate(wk(2*mxl-1,2*mxl-1),stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'uutrns','allocation fails')
      if(il .ne. 1 .and. il .ne. 2)
     &                call errtrp(1,'uutrns','illegal il')
      ir=3-il
c     if mxl=1, there is no need to transform anything.
      do 10 l=2,mxl
      ms=2*l-1
      mi=mxl-l
      do 10 iop=1,24
c     urotat(1,1,iop)=0 means the operation corresponding to iop
c     does not exist, e.g. for the hexagonal case iop runs only
c     from 1 to 12.
      if(abs(urotat(1,1,iop)) .gt. 1d-10) then
      call clrarc(wk,(2*mxl-1)**2)
      do 20 mc=1,ms 
      mcu=mi+mc
      do 20 mr=1,ms
      mru=mi+mr
      if(abs(u(mru,mcu,ir)) .gt. zero) then
      mb=ms*(mr-1)
      do 30 m=1,ms
      mmr=mb+m
   30 wk(m,mc)=wk(m,mc)+urotat(mmr,l,iop)*u(mru,mcu,ir)
      endif
   20 continue
      call clrarc(urotat(1,l,iop),ms**2)
      do 40 mc=1,ms
      mcu=mi+mc
      do 40 mr=1,ms
      mru=mi+mr
      if(abs(u(mru,mcu,il)) .gt. zero) then
      do 50 m=1,ms
      mrm=ms*(m-1)+mr
   50 urotat(mrm,l,iop)=urotat(mrm,l,iop)+u(mru,mcu,il)*wk(mc,m)
      endif
   40 continue
      endif
   10 continue
      deallocate(wk,stat=ierr)
      end
