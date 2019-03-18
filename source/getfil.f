
      subroutine getfil(lunit,filnam)
c----------------------------------------------------------------------
c     This program opens file with unit=lunit and file=filnam.
c     When file does not exist, a new file is created and eof
c     is marked.
c     coded by H.Akai, 1988, Tokyo(ISSP)
c     revised by H.Akai, Feb. 1993, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character filnam*(*)
      logical exs,opn
      inquire(file=filnam,exist=exs,opened=opn,number=ifil)
      if(opn) close(ifil)
      open(lunit,file=filnam,form='unformatted',status='unknown')
      write(*,'(2a)') '   file to be accessed=',filnam
      if(exs) return
c     if the eof is not properly detected, activate the following
c     two statements.
c     endfile(lunit)
c     rewind(lunit)
      write(*,'(a)') '   created'
      return
      end
