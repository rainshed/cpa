      subroutine errtrp(idsp,routin,msg)
c-----------------------------------------------------------------------
c     Error trap:
c     idisp=0...error causing no imediate terminataion
c     idisp=1...error causing imediate termination
c     idisp=2...give warning message
c     idisp=3...give some usuful information
c     Coded by H.Akai, 16 Oct. 1993, Osaka
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      character*(*) routin,msg
      go to (10,20,30,40), idsp+1
      write(*,'(a)')' ***err in errtrp...unexpected idsp'
      stop
   10 write(*,'(4a)')' ***err in ',routin,'...',msg
      return
   20 write(*,'(4a)')' ***err in ',routin,'...',msg
      stop
   30 write(*,'(4a)')' ***wrn in ',routin,'...',msg
      return
   40 write(*,'(4a)')' ***msg in ',routin,'...',msg
      return
      end
