      subroutine wrtspc(spctrl,vkp,is,e,mse,kcrt,kblst,nk1,nk3,ef
     &                 ,unit,iwrt)
c-----------------------------------------------------------------------
c     Given array data {{spctrl(k,kp), k=1,mse},kp=1,nk3} and
c     {e(k),k=1,mse}, this program write down plot data for the
c     dispersion relation E(k). The array {kcrt(k),k=1,kblst}
c     contains the pointer for the k-points that are on a symmetry
c     point. The data are out put on a file unit=27 and 28 for spin
c     up and down, respectively.
c
c     Coded by H. Akai, 6 June 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 spctrl(mse,nk3),e(mse,2)
      real*8 vkp(3,nk1+nk3),ef(2)
      integer kcrt(kblst)
      data small/1d-5/
      if(iwrt .eq. 1) then
c     --- output with format (a)
      sftfct=5d-1*unit
      estep=dble(e(2,is))-dble(e(1,is))
c     write(*,'(1x,a,i2)')'#  A(E,k) for spin =',is
      write(26+is,'(a,i2)')'#  A(E,k) for spin =',is
      do kp=1,nk3
         kk=nk1+kp
c     write(*,'(/1x,a,3f12.5)')'   k=',(vkp(i,kk),i=1,3)
c     write(27,'(/a,3f12.5)')'#  k=',(vkp(i,kk),i=1,3)
         if(kp .eq. 1) then
            dist=0d0
         else
            dist=dist+sqrt( (vkp(1,kk)-vkp(1,kk-1))**2
     &           +(vkp(2,kk)-vkp(2,kk-1))**2
     &           +(vkp(3,kk)-vkp(3,kk-1))**2 )*sftfct
         endif
c--- output in (x,y,z) form ----
         write(26+is,*)
         do k=2,mse
            write(26+is,'(1x,3f15.7)')dist,5d-1*dble(e(k,is)+e(k-1,is))
     &           -ef(is),max(0d0,dimag((spctrl(k,kp)
     &           -spctrl(k-1,kp))/(e(k,is)-e(k-1,is))))
         end do 
      end do
      return
      else if(iwrt .eq. 2) then
c     --- output with format (b)
c--- output in matrix form (for Igor)----
c    horizontal: energy, vertical: k-point
      do 20 kp=1,nk3
   20 write(26+is,'(1x,500f14.7)')(max(0d0,dimag((spctrl(k,kp)
     &    -spctrl(k-1,kp))/(e(k,is)-e(k-1,is)))),k=2,mse)
      return
      else if(iwrt .eq. 3) then
c     --- output with format (c)
c--- output in matrix form (for gnuplot)----
c    horizontal: k-point, vertical: energy
      nkd=max(3,nk3/60)
      do 30 k=1,mse-1
      do 40 kp=1,nk3
      spctrl(k,kp)=max(0d0,dimag((spctrl(k+1,kp)
     &              -spctrl(k,kp))/(e(k+1,is)-e(k,is))))
      do 40 kc=2,kblst-1
   40 if(kp .eq. kcrt(kc)) spctrl(k,kp)=250d0
      if(abs(dble(e(k,is))-ef(is)) .lt. small) then
      do 50 kp=1,nk3,nkd
   50 spctrl(k,kp)=250d0
      endif
   30 write(26+is,'(1x,500f14.7)')(dble(spctrl(k,kp)),kp=1,nk3)
      return
      else
      call errtrp(1,'wrtspc','illegal iwrt')
      endif
      end
