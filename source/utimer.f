      subroutine utimer(time,i)
c-----------------------------------------------------------------------
c     Driving timer routine.
c     coded by H.Akai, 1980, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      time1=time
      call uclock(time)
      if(i .eq. 0) return
      write(6,1100)i,time-time1
 1100 format('   interval=',i3,'     cpu time=',f10.2,' sec')
      return
      end
