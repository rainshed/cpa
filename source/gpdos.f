      program dosplt
c-----------------------------------------------------------------------
c     construct DOS data from cpa98 output.
c     coded by H.Akai, 22 Feb. 1998, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(msex=201,icmpx=10,maxlin=10000,mmxl=4)
      real*8 dat(mmxl*icmpx+3,msex)
      character file*256,a*120,fmt*80,fmtx*80
      logical exist
c--- get input file
      file=''
      call getarg(1,file)
   50 if(file .eq. '') then
      write(*,'(a)')' input file name ?'
      read(*,'(a)',end=10)file
      endif
      call lftfil(file)
      call chleng(file,len)
      if(len .le. 0) go to 10
c
c     If a directory of the same name exists, inquire can not
c     ditect existence of a file correctly. Therefore one should
c     use open statement with error return.
      inquire(file=file,exist=exist)
      if(.not. exist) then
      write(*,'(a)')' file not found'
      file=''
      go to 50
      endif
      open(12,file=file)
c---  create output file
      ip=index(file,'.')
      if(ip .ge. 2) file=file(1:ip-1)
      call lftfil(file)
      call chleng(file,len)
      if(len .le. 0) go to 10
c
c--- read-in files (unit=12) and find DOS data.
c
      mse=0
      mxl=0
      icmp=0
      do 70 line=1,maxlin
c     --- get the number of meshes, l_max, and number of components
      read(12,'(a)',end=10)a
      if(a(5:28) .eq. 'meshr  mse    ng     mxl') then
      read (12,*)meshr,mse,ng,mxl
      endif
      if(a(4:8) .eq. 'ntyp='.and. a(13:17) .eq. 'natm='
     &   .and. a(22:27) .eq. 'ncmpx=') then
      read(a(28:29),'(i2)')icmp
      endif
c     --- the last data before the DOS data begin are adopted
      if(a(2:19) .eq. 'DOS of component 1') go to 100
   70 continue
      write(*,'(a)')' seems not to contain DOS data'
      go to 10
  100 write(*,'(a,3i3)')' mse,mxl,icmp=',mse,mxl,icmp
c     --- check the sizes ---
      if(mse .gt. msex) then
      write(*,'(a,i3)')' ***error...mse>',msex
      go to 10
      endif
      if(mxl .gt. mmxl) then
      write(*,'(a,i2)')' ***error...mxl>',mmxl
      go to 10
      endif
      if(icmp .gt. icmpx) then
      write(*,'(a,i2)')' ***error...icmp>',icmpx
      go to 10
      endif
c     --- Now the format statement can be constructed
      write(fmtx,'(a,i1,a)')'(1x,f7.4,3x,',mxl,'f10.4)'
c     write(*,*)'format=',fmtx
      open(13,file=file(1:len)//'.plt',status='unknown')
c     write(13,'(a)')'#!/usr/bin/gnuplot -persist'
      write(13,'(a)')'#!/usr/bin/gnuplot'
      write(13,'(2a)')'# set terminal postscript landscape noenhanced ',
     &               'monochrome dashed defaultplex "Helvetica" 14'
      write(13,'(a)')'# set output "dos.eps"'
      write(13,'(a)')'set term x11'
      write(13,'(a)')'set border 15 lw 0.1'
      write(13,'(a)')'set size ratio 0.7071'
      write(13,'(a)')'set xzeroaxis lt -1 lw 0.1'
      write(13,'(a)')'set yzeroaxis lt -1 lw 0.1'
      write(13,'(a)')'set mxtics 2'
      write(13,'(a)')'set mytics 2'
      write(13,'(a)')'set xtics border mirror norotate 0.2'
      write(13,'(a)')'set ytics border mirror norotate 10'
      write(13,'(3a)')'set title "',file(1:len),
     &                 '" offset 0.000000,0.000000'
      write(13,'(2a)')
     & 'set xlabel "Energy relative to Fermi energy (Ry)" ',
     & ' offset 0.000000,0.000000'
      write(13,'(a)')
     & 'set xrange [ * : * ] noreverse nowriteback'
      write(13,'(a)')
     & 'set ylabel "DOS (1/Ry/unit cell)" offset 0.000000,0.000000'
      write(13,'(a)')
     & 'set yrange [ * : * ] noreverse nowriteback'
      write(13,'(2a)')'plot "-" u 1:2 t "" w l ls 1, ',
     &             '"-" u 1:(-1*$2) t "" w l ls 1'
      do 20 is=1,2
c	write(*,'(1x,i3)')is
c	write(*,'(1x,a)')file(1:len)//ext(is)
      do 30 line=1,maxlin
      if(a(2:19) .eq. 'DOS of component 1') then
c     write(*,*)'DOS           1 met'
      read(12,fmtx)((dat(l,k),l=3,mxl+3),k=1,mse)
c     write(*,fmtx)((dat(l,k),l=3,mxl+3),k=1,mse) 
   60 continue
      do 40 i=2,icmp
c     write(*,*)'DOS',i,' met'
      i0=mxl*(i-1)+4
      read(12,'(//)')
   40 read(12,fmtx)
     &    (dmy,(dat(l,k),l=i0,i0+mxl-1),k=1,mse)
      endif
      if(a(2:10) .eq. 'total DOS') then
c     write(*,*)'total DOS met'
      read(12,'(1x,f12.7,3x,f13.5)')(dat(1,k),dat(2,k),k=1,mse-1)
      dat(1,mse)=dat(1,mse-1)
      dat(2,mse)=dat(2,mse-1)
      write(fmt,'(a,i2,a)')'((f8.4,f7.2,f8.4,',mxl*icmp,'f7.2))'
      write(13,fmt)((dat(l,k),l=1,mxl*icmp+3),k=1,mse)
      write(13,'(a)')'end'
      go to 20
      endif
   30 read(12,'(a)',end=10)a
   20 read(12,'(a)',end=10)a
      write(13,'(a)')'pause -1'
c     write(13,'(a)')'#    EOF'
      close(13)
      call system('gnuplot '//file(1:len)//'.plt')
   10 continue
      end
      subroutine lftfil(fil)
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character fil*(*)
      n=len(fil)
      j=0
      is=1
      do 30 i=1,n
      is=i
      if(llt(fil(i:i),'0') .or. llt('9',fil(i:i))) go to 40
   30 continue
   40 do 10 i=is,n
      if(fil(i:i).eq.' ') go to 10
      j=j+1
      fil(j:j)=fil(i:i)
   10 continue
      if(j.ge.n) return
      fil(j+1:n)=' '
      return
      end
      subroutine chleng(a,ln)
c---------------------------------------------------------------------
      character a*(*)
      n=len(a)
      do 10 i=n,1,-1
      ln=i
      if(a(i:i) .ne. ' '.and. ichar(a(i:i)) .gt. 27) return
   10 continue
      ln=0
      return
      end
