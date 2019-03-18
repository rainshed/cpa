      subroutine udate(today)
c-----------------------------------------------------------------------
c     This program returns the character date corresponding
c     the today's date.
c     This is a machin dependent program.
c     coded by H. Akai 3 Nov 1996, Osaka
c-----------------------------------------------------------------------
      character today*11,mname(12)*4,day*8,time*10,zone*5
      data mname/'Jan','Feb','Mar','Apr','May','Jun','Jul',
     &           'Aug','Sep','Oct','Nov','Dec'/ 
c
c     --- for hp ---
c     call date(today)
c     return
c
c     --- for digital unix and ibm RISC 6000 ---
      call date_and_time(day,time,zone)
      read(day,'(i4,i2,i2)')iyear,month,iday
      write(today,'(i2,a,a3,a,i4)')
     &   iday,'-',mname(month),'-',iyear
c
c     --- for MS-FOTRAN ---
c     call getdat(iyear,month,iday)
c     write(today,'(i2,a,a3,a,i4)')
c    &   iday,'-',mname(month),'-',iyear
c      return
      end
