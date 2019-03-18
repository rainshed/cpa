      subroutine rotave(g,urotat,isymop,kmx,mxl,lmxtyp,ncmpx)
c-----------------------------------------------------------------------
c     Rotate the matrix f(mr,mc,iatm) for each atom according to the
c     symmetry table given by urotat, irotat and isymop, and add up
c     the resulting matrix elements.
c
c     urotat: rotation matrices for any 24 cubic or 12 hexgonal
c             operations (inversion is omitted).
c     isymop: if isymop(iop)=1, then rotation specified by iop is
c             allowed from the symmetry.
c
c     Coded by H.Akai, 28 Nov. 2014
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(ismx=1)
      complex*16 g(ismx*mxl**2,ncmpx,kmx),urotat((2*mxl-1)**2,mxl,24),u
      complex*16,allocatable::wk(:,:)
      integer isymop(24,2),lmxtyp(*)
      data zero/1d-10/
      allocate(wk(mxl**2,ncmpx))
      do 10 is=0,ismx-1
      do 10 k=1,kmx
      call clrarc(wk,mxl**2*ncmpx)
      weight=0d0
      do 20 iop=1,24
c     --- only rotations compatible with the crystal symmetry
c         are considered.
      if(isymop(iop,1) .eq. 1 .or. isymop(iop,2) .eq. 1) then
c     Since no difference between proper or inproper rotation arizes for
c     the rotation of a diagonal matrix, the proper ratation is
c     applied irrespective of the real symmetry.
      weight=weight+1d0
      do 30 l=1,mxl
      mx=2*l-1
      lb=(l-1)**2
      do 30 m1=1,mx
      do 30 m2=1,mx
      m1m2=mx*(m2-1)+m1
      if(abs(urotat(m1m2,l,iop)) .gt. zero) then
c     --- only non-zero elements of the rotation matrix are
c         taken into account.
      u=urotat(m1m2,l,iop)*conjg(urotat(m1m2,l,iop))
      do 40 i=1,ncmpx
      call jipinv(i,ityp,icmp)
      if(l .le. lmxtyp(ityp)+1) then
      mmx=(lmxtyp(ityp)+1)**2
      wk(lb+m1,i)=wk(lb+m1,i)+g(lb+m2+mmx*is,i,k)*u
      endif
   40 continue
      endif
   30 continue
      endif
   20 continue
      factor=1d0/weight
c     if(k .eq. kmx) then
c     write(*,*)'---------------'
c     write(*,*)'weight=',weight
c     write(*,*)'f'
c     write(*,'((1x,1p9e20.13))')(dble(g(m,1,k)),m=1,9)
c     write(*,*)'wk'
c     write(*,'((1x,1p9e20.13))')(factor*dble(wk(m,1)),m=1,9)
c     write(*,*)
c     endif
      do 10 i=1,ncmpx
      call jipinv(i,ityp,icmp)
      mmx=(lmxtyp(ityp)+1)**2
      do 10 m=1,mmx
   10 g(m+mmx*is,i,k)=factor*wk(m,i)
      deallocate(wk)
      end
