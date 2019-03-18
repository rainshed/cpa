      subroutine setomp
c-----------------------------------------------------------------------
c     get the omp wall time and set number of threads
c     coded by H. Akai 10 Aug. 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      save tstart,nthread
      tstart=omp_get_wtime()
      nthread=numcor()
      if(nthread .gt. 16) nthread=nthread-1 
c     nthread=max(numcor()-1,1)
      call omp_set_num_threads(nthread)
c     call omp_get_num_threads(nthread)
      return
      entry endomp
c-----------------------------------------------------------------------
c     get the end wall time and output elapsed time.
c-----------------------------------------------------------------------
      tend=omp_get_wtime()
      write(*,'(1x,a,t15,f12.2,a,i3,a//)')
     & 'elapsed time ',tend-tstart,' sec (',nthread,' threads)'
      end
