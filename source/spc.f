      program spc
c-----------------------------------------------------------------------
c     construct gnuplot batch file and submit it.
c     coded by H.Akai, 13 June 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character file*256,tmpfil*6,command*256
      logical exist
c--- get input file
      file=' '
      call getarg(1,file)
      if(file .eq. ' ') then
      write(*,'(a)')'Usage: spc file'
      stop
      endif
      len=lftfil(file)
c
c     If a directory of the same name exists, inquire cannot
c     ditect existence of a file correctly. Therefore one should
c     use 'open' statement together with error return.
      inquire(file=file,exist=exist)
      if(.not. exist) then
      write(*,'(a/a)')'File not found','Usage: spc file'
      stop
      endif
c---  create a temporary file
      do 10 i=0,999
      write(tmpfil,'(a,i3.3)')'tmp',i
      inquire(file=tmpfil,exist=exist)
      if(.not. exist) go to 20
   10 continue
      write(*,'(a)')'create temporary file failed'
      go to 30
   20 open(13,file=tmpfil,status='new')
c     --- construct a gnuplot program
c     write(13,'(a)')'#!/usr/bin/gnuplot -persist'
      write(13,'(a)')'#!/usr/bin/gnuplot'
      write(13,'(2a)')'set terminal postscript landscape noenhanced ',
     &               'monochrome dashed defaultplex "Helvetica" 14'
      write(13,'(a)')'#set output "spc.eps"'
      write(13,'(a)')'set term x11'
      write(13,'(a)')'set pm3d map'
      write(13,'(a)')'unset tics'
      write(13,'(a)')
     & 'set palette defined (0 "white", 60 "blue", 259 "blue")'
      write(13,'(3a)')'splot "',file(1:len),'" matrix'
c     write(13,'(a)')'#pause -1'
      write(13,'(a)')'pause -1'
      call system('gnuplot '//tmpfil)
   30 close(13)
      call system('rm '//tmpfil)
      end
      function lftfil(fil)
c---------------------------------------------------------------------
c     Get rid of the blank from the string. Remaining part is
c     compressed to the left. The function returns the length of
c     the non-zero part.
c     coded by H.Akai, 1986, Juelich
c     Modified by H.Akai, 2007
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character fil*(*)
      n=len(fil)
      j=0
      do 10 i=1,n
      if(fil(i:i) .eq. ' ') go to 10
      j=j+1
      fil(j:j)=fil(i:i)
   10 continue
      lftfil=j
      if(j .ge. n) return
      fil(j+1:n)=' '
      return
      end
      subroutine chleng(a,ln)
c---------------------------------------------------------------------
c     Returns the length of the non blank part of the character a.
c     coded by H.Akai
c---------------------------------------------------------------------
      character a*(*)
      n=len(a)
      do 10 i=n,1,-1
      ln=i
      if(a(i:i) .ne. ' ') return
   10 continue
      ln=0
      return
      end
