      subroutine bzmrot(f,natm,irotat,urotat,isymop,kmx
     &           ,convrg,mxl,snor,lmxtyp,itype,ls)
c-----------------------------------------------------------------------
c     Rotate the matrix f(mr,mc,iatm) for each atom according to the
c     symmetry table given by urotat, irotat and isymop, and add up
c     the resulting matrix elements.
c
c     urotat: rotation matrices for any 24 cubic or 12 hexgonal
c             operations (inversion is omitted).
c     irotat: associated rotation in the atomic position. If atom 2
c             comes to the position where atom 1 was sitting before the
c             rotation, then irotat(1)=2, etc.
c     isymop: if isymop(iop)=1, then rotation specified by iop is
c             allowed from the symmetry.
c
c     Coded by H.Akai, 29 July 1997, Muenchen
c     Modified by H. Akai, 20 Aug. 2009
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(ismx=1)
      complex*16 f(ismx*mxl**2,ismx*mxl**2,natm,kmx)
     &          ,u,urotat((2*mxl-1)**2,mxl,24)
      complex*16,allocatable::wk1(:,:,:),wk2(:,:,:)
      integer irotat(natm,24,2),isymop(24,2),itype(natm),lmxtyp(*)
      logical convrg(kmx)
      data zero/1d-10/
      allocate(wk1(mxl**2,mxl**2,natm),wk2(mxl**2,mxl**2,natm))
c     return
c     do 100 iop=1,24
c     write(*,'(3i2)')iop,isymop(iop,1),isymop(iop,2)
c     if(isymop(iop,1) .eq. 1) then
c     write(*,'(4i3)')(irotat(i,iop,1),i=1,natm)
c     elseif(isymop(iop,2) .eq. 1) then
c     write(*,'(4i3)')(irotat(i,iop,2),i=1,natm)
c     endif
c 100 continue
c     diff=0d0
      do 10 is1=0,ismx-1
      do 10 is2=0,ismx-1
      do 10 k=1,kmx
c     --- calculation is needed only when we have not yet gotten
c         a converged results.
      if(.not. convrg(k)) then
      call clrarc(wk1,mxl**4*natm)
      weight=0d0
      do 20 ip=1,2
      do 20 iop=1,24
c     --- only rotations compatible with the crystal symmetry
c         are considered.
      if(isymop(iop,ip) .eq. 1) then
c     write(*,'(1x,2i3,2x,20i3)')iop,isymop(iop),
c    &             (irotat(i,iop),i=1,natm)
      weight=weight+1d0
      call clrarc(wk2,mxl**4*natm)
c     --- we first calculate wk2=f*U^+, where U^+ is the Hermite
c         conjugate of the rotation matrix U.
      do 30 l=1,mxl
      p=(-1d0)**((ip-1)*(l-1))
c     ll=l
c     call dspu(urotat(1,l,iop),ll)
      mx=2*l-1
      lb=(l-1)**2
      do 30 m1=1,mx
      m0=mx*(m1-1)
      do 30 m2=1,mx
      m2m1=m0+m2
      if(abs(urotat(m2m1,l,iop)) .gt. zero) then
c     --- only non-zero elements of the rotation matrix are
c         taken into account. In order to handle this procedure
c         efficiently, the order of summation appearing in
c         the matrix product is changed from a straight forward
c         way.
c     --- take account of the effects of the inversion.
      u=conjg(p*urotat(m2m1,l,iop))
      mr=lb+m1
      mc=lb+m2
      do 40 i=1,natm
c     j=irotat(i,iop,ip)
      if(l .le. lmxtyp(itype(i))+1) then
      mmx=(lmxtyp(itype(i))+1)**2
      do 42 m=1,mmx
   42 wk2(m,mc,i)=wk2(m,mc,i)+f(m+mmx*is1,mr+mmx*is2,i,k)*u
      endif
   40 continue
      endif
   30 continue
c     --- next, we calculate U*wk2 and the results are accumulated
c         on wk1.
      do 50 l=1,mxl
      p=(-1d0)**((ip-1)*(l-1))
      mx=2*l-1
      lb=(l-1)**2
      do 50 m1=1,mx
      m0=mx*(m1-1)
      do 50 m2=1,mx
      m2m1=m0+m2
      if(abs(urotat(m2m1,l,iop)) .gt. zero) then
      u=p*urotat(m2m1,l,iop)
      mc=lb+m1
      mr=lb+m2
      do 60 i=1,natm
c     --- we also have to take it account that the atom positions
c         are rotated.
      if(l .le. lmxtyp(itype(i))+1) then
      mmx=(lmxtyp(itype(i))+1)**2
      j=irotat(i,iop,ip)
      do 62 m=1,mmx
   62 wk1(mr,m,i)=wk1(mr,m,i)+u*wk2(mc,m,j)
      endif
   60 continue
      endif
   50 continue
      endif
   20 continue
c     --- now the allowed rotations have been taken, and we
c         can replace the old f by a new one. Also the summation
c         with respect to -k, which was omitted from the BZ mesh
c         for the cases of
c           i) time reversal symmetry exists
c           ii) inversion symmetry exists
c         are now taken into account.
c         "snor" is a presumed normalization factor.
c     go to 100
      if(ls .eq. 0) then
      weight=2d0*weight
c     --- case i)
      call equarc(wk1,wk2,mxl**4*natm)
      do 70 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      do 70 mc=1,mmx
      do 70 mr=1,mmx
   70 wk1(mr,mc,i)=wk2(mr,mc,i)+wk2(mc,mr,i)
      else if(isymop(1,2) .eq. 1) then
c     --- case ii)
      weight=2d0*weight
      call equarc(wk1,wk2,mxl**4*natm)
      do 80 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      j=irotat(i,1,2)
      do 80 mc=1,mmx
      l1=lindx(mc)-1
      do 80 mr=1,mmx
      l2=lindx(mr)-1
   80 wk1(mr,mc,j)=wk2(mr,mc,j)+wk2(mr,mc,i)*(-1d0)**(l1+l2)
      endif
c 100 continue
c     --- do nothing if neither case i) nor case ii)
      factor=1d0/weight/snor
c     if(k .eq. kmx) then
c     write(*,*)'---------------'
c     write(*,*)'f_1'
c     write(*,'((1x,1p9e20.13))')
c    &   ((dble(f(mr,mc,1,k)),mc=1,9),mr=1,9)
c     write(*,*)'f_2'
c     write(*,'((1x,1p9e20.13))')
c    &   ((dble(f(mr,mc,2,k)),mc=1,9),mr=1,9)
c     write(*,*)'wk1_1'
c     write(*,'((1x,1p9e20.13))')
c    &   ((dble(factor*wk1(mr,mc,1)),mc=1,9),mr=1,9)
c     write(*,*)'wk1_2'
c     write(*,'((1x,1p9e20.13))')
c    &   ((dble(factor*wk1(mr,mc,2)),mc=1,9),mr=1,9)
c     write(*,*)
c     endif
      do 90 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      do 90 mc=1,mmx
      do 90 mr=1,mmx
c     diffx=abs(f(mr+mmx*is1,mc+mmx*is2,i,k)-factor*wk1(mr,mc,i))
c     diff=max(diff,diffx)
c     if(diffx .gt. 1d1) then
c     write(*,*)k,i,mc,mr
c     write(*,'(1p2e20.13)')f(mr+mmx*is1,mc+mmx*is2,i,k)
c     write(*,'(1p2e20.13)')factor*wk1(mr,mc,i)
c     endif
   90 f(mr+mmx*is1,mc+mmx*is2,i,k)=factor*wk1(mr,mc,i)
c  90 continue
      endif
   10 continue
c     write(*,*)'diff=',diff
c     do 1010 l=1,natm
c     write(*,*)'atom=',l
c1010 write(*,'(1x,6f12.7)')((f(i,j,l,kmx),j=2,4),i=2,4)
      deallocate(wk1,wk2)
      end
