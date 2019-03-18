      function numcor()
c-----------------------------------------------------------------------
c     This function returns the number of cores of the system.
c     Coded by H. Akai, 13 Feb. 2015, Tokyo
c     Modified to be also adapted for MAC OSX
c     by H. Akai 8 Aug. 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mxlin=10000)
      character*80 file,a,tmpfil*6,command
      logical exist
      data inicor/4/
c     ---open /proc/cpuinfo to get processor information
      file='/proc/cpuinfo'
      inquire(file=file,exist=exist)
      if(exist) then
      open(13,file=file,action='read')
c     ---read in /proc/cpuinfo (unit=13)
      numcor=0
      do 10 line=1,mxlin
      read(13,'(a80)',end=20)a
   10 if(index(a,'processor') .gt. 0)numcor=numcor+1
      call errtrp(1,'numcor','mxlin too small')
   20 close(13)
      else
c     ---try if the sysmtem is Mac OSX
      do 30 i=0,999
      write(tmpfil,'(a,i3.3)')'tmp',i
      inquire(file=tmpfil,exist=exist)
      if(.not. exist) go to 40
   30 continue
      write(*,'(a)')' create temporary file failed'
      go to 50
   40 open(13,file=tmpfil)
c     ---the following command is available in Mac OSX
      command='sysctl -n machdep.cpu.thread_count>'//tmpfil
      istat=system(command)
      if(istat .eq. 0) then
      read(13,'(i3)',end=50,err=50)numcor
      close(13)
      command='rm '//tmpfil
      istat=system(command)
      return
      endif
   50 continue
      write(*,'(a,i2,a)')' cpuinfo not found: number of cores='
     &     ,inicor,' is assumed.'
      numcor=inicor
      return
      endif
      end
