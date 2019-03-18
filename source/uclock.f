      subroutine uclock(time)
c-----------------------------------------------------------------------
c     User clock which initializes clock when called first time
c     and gives lap time after the first call. The program is
c     maching dependent.
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical first
      character today*11,date*8,now*10,zone*5
      integer ifval(8)
      save time0,jday0,first
      data time0/0d0/, first/.true./
c
c     --- for ibm vs ---
c     real*4 ts
c     call frist(ts)
c     t=-dble(ts)*1d-2
c     if(first) then
c     time0=t
c     first=.false.
c     endif
c     time=t-time0
c     return
c
c     --- for vpp/500  ---
c     call clockm(ntime)
c     t=dble(ntime)*1d-3
c     if(first) then
c     time0=t
c     first=.false.
c     endif
c     time=t-time0
c     return
c     call clockm(ntime)
c     t=dble(ntime)*1d-3
c     if(first) then
c     time0=t
c     first=.false.
c     endif
c     time=t-time0
c     return
c
c     --- for nec acos ---
c     call cptime(ntime)
c     t=dble(ntime)*1d-3
c     if(first) then
c     time0=t
c     first=.false.
c     endif
c     time=t-time0
c     return
c
c     --- for nec sx ---
c     call clock(t)
c     if(first) then
c     time0=t
c     first=.false.
c     endif
c     time=t-time0
c     return
c
c     --- for cray ---
c     call second(t)
c     if(first) then
c     time0=t
c     first=.false.
c     endif
c     time=t-time0
c     return
c
c     --- for hitac s820 ---
c     if(first) then
c     call xclock
c     first=.false.
c     t=time0
c     return
c     endif
c     call xclock(t,5)
c     time=t-time0
c     return
c
c     --- for hp and old dec ---
c     t=dble(secnds(0e0))
c     call idate(month,iday,iyear)
c     jday=julday(month,iday,iyear)
c     if(first) then
c     call udate(today)
c     write(*,'(1x,a)')today
c     time0=t
c     jday0=jday
c     first=.false.
c     endif
c     time=t-time0+86400d0*dble(jday-jday0)
c     return
c
c      --- for digital unix and ibm RISC 6000 ---
c      t=dble(secnds(0e0))
      call cpu_time(t)
      call date_and_time(date,now,zone)
      read(date,'(i4,i2,i2)')iyear,month,iday
      read(now,'(i2,i2,i2)')ihour,minute,isec
      jday=julday(month,iday,iyear)
      if(first) then
      call udate(today)
      write(*,'(1x,a/4x,i2,a,i2.2,a,i2.2)')
     &      today,ihour,':',minute,':',isec
      time0=t
      jday0=jday
      first=.false.
      endif
      time=t-time0+86400d0*dble(jday-jday0)
      return
c
c     --- for MS-FOTRAN ---
c     call gettim(ihr,imin,isec,i100th)
c     t=(dble(ihr)*60d0+dble(imin))*60d0+dble(i100th)*1d-2
c     call getdat(iyear,month,iday)
c     jday=julday(month,iday,iyear)
c     if(first) then
c     call date(today)
c     write(*,'(1x,a)')today
c     time0=t
c     jday0=jday
c     first=.false.
c     endif
c     time=t-time0+86400d0*dble(jday-jday0)
c     return
      end
