      subroutine sbtime(isub,icount)
c-----------------------------------------------------------------------
c     Count cpu time used in a specific subroutine.
c     To initialize the couter, use
c          call sbtime(0,0)
c     Without the initialization, any call to this routine simply causes
c     an immediate return.
c     To accumulate the cputime used in a subroutine to the counter, use
c          call sbtime(isub,0)
c     at the bigining of the subroutine, where isub distinguishes
c     the subroutine (any unique positibve integer number less than
c     isubmx in the parameter list.
c     At the end of the same routine (just before the return), use
c           call sbtime(isub,n)
c     where n can be either 1 or any integer number that will be added
c     to the number counter ncount.
c     The results of counting is displayed by
c           call sbtime(-1,0)
c     Here, -1 could be any negative integer and the second argument
c     cound be any number. The counter is inactivated at the same
c     time, sleeping until it is again activated by initialization.
c
c     coded by H. Akai, 31. Mar. 02, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(isubmx=20)
      real*8 tcount(isubmx),tmp(isubmx)
      integer ncount(isubmx)
      logical skip
      save tcount,tmp,ncount,skip
      data skip/.true./
      if(isub .eq. 0) then
      skip=.false.
      do 10 i=1,isubmx
      tcount(i)=0d0
   10 ncount(i)=0
      return
      endif
      if(skip) return
      if(isub .gt. isubmx) call errtrp(1,'sbtime','isub too large')
      if (isub .gt. 0) then
      call cpu_time(time)
      if(icount .eq. 0) then
      tmp(isub)=time
      return
      else
      tcount(isub)=tcount(isub)+time-tmp(isub)
      ncount(isub)=ncount(isub)+icount
      return
      endif
      else
      im=0
      do 20 i=1,isubmx
   20 if(ncount(i) .gt. 0) im=i
      irp=im/10
      ied=mod(im,10)
      write(*,*)
      write(*,'(1x,a)')'sbtime report'
      do 30 ic=1,irp
      ib=10*(ic-1)
      write(*,'(1x,a10,10i10)')'routine  ',(i,i=ib+1,ib+10)
      write(*,'(1x,a10,10i10)')'count    ',(ncount(i),i=ib+1,ib+10)
      write(*,'(1x,a10,10f10.2)')'cpu(sec) ',(tcount(i),i=ib+1,ib+10)
   30 write(*,*)
      ib=10*irp
      write(*,'(1x,a10,10i10)')'routine  ',(i,i=ib+1,im)
      write(*,'(1x,a10,10i10)')'count    ',(ncount(i),i=ib+1,im)
      write(*,'(1x,a10,10f10.2)')'cpu(sec) ',(tcount(i),i=ib+1,im)
      write(*,*)
      skip=.true.
      return
      endif
      end
