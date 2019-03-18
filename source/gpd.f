      program gpd
c-----------------------------------------------------------------------
c     construct DOS data from cpa98 output.
c     coded by H.Akai, 22 Feb. 1998, Osaka
c     modifed by H. Akai, 27 Nov. 2010
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(maxlin=10000)
      real*8,allocatable::dat(:,:,:)
      character file*256,a*120,fmt*80,fmtx*80,magtyp*4
      logical exist
      data big/1d30/
c--- get input file
      file=' '
      call getarg(1,file)
      if(file .eq. ' ') then
      write(*,'(a)')'Usage: gpd file'
      stop
      endif
      call lftfil(file)
c
c     If a directory of the same name exists, inquire cannot
c     ditect existence of a file correctly. Therefore one should
c     use 'open' statement together with error return.
      inquire(file=file,exist=exist)
      if(.not. exist) then
      write(*,'(a/a)')'File not found','Usage: gpd file'
      stop
      endif
      open(12,file=file)
c---  create output file
      ip=index(file,'.')
      if(ip .ge. 2) file=file(1:ip-1)
      call lftfil(file)
      call chleng(file,len)
      if(len .le. 0) then
      write(*,'(a/a)')'Illegal file name','Usage: gpd file'
      close(12)
      stop
      endif
c
c--- read-in files (unit=12) and find DOS data.
c
      mse=0
      mxl=0
      icmp=0
      magtyp=' '
      do 50 line=1,maxlin
c     --- get the number of meshes, l_max, and number of components
      read(12,'(a)',end=10)a
      if(index(a,'meshr')*index(a,'mse')*index(a,'ng')*index(a,'mxl')
     &    .gt. 0) read (12,*)meshr,mse,ng,mxl
      n=index(a,'ncmpx=')
      if(n .gt. 0) read(a(n+6:n+9),'(i3)')icmp
c     --- get magtyp
      n=index(a,'magtyp=')
      if(n .gt. 0) magtyp=a(n+7:n+10)
c     --- the last data before the DOS data begin are adopted
      if(index(a,'DOS of component 1') .gt. 0) go to 100
   50 continue
      write(*,'(a)')'Seems not to contain DOS data.'
      close(12)
      stop
c 100 write(*,'(a,3i3)')' mse,mxl,icmp=',mse,mxl,icmp
  100 continue
      if(magtyp .eq. 'mag') then
      ismx=2
      else if(magtyp .eq. 'nmag') then
      ismx=1
      else
      write(*,'(a,a4)')' ***error...illegal magtyp=',magtyp
      close(12)
      stop
      endif
      allocate(dat(mxl*icmp+3,mse,2),stat=ierr)
      if(ierr .ne. 0) then
      write(*,'(a)')' ***err in gpd...allocation fails'
      stop
      endif
c     write(*,*)'ismx=',ismx
c     --- Now the format statement can be constructed
      write(fmtx,'(a,i1,a)')'(1x,f7.4,3x,',mxl,'f10.4)'
      write(fmt,'(a,i2,a)')'((f8.4,f7.2,f8.4,',mxl*icmp,'f7.2))'
c     write(*,*)'format=',fmtx
c     --- read in DOS data
      do 20 is=1,ismx
      do 30 line=1,maxlin
      if(index(a,'DOS of component 1') .gt. 0) then
c     write(*,*)'DOS           1 met'
      read(12,fmtx)((dat(l,k,is),l=3,mxl+3),k=1,mse)
c     write(*,fmtx)((dat(l,k),l=3,mxl+3,is),k=1,mse) 
      do 40 i=2,icmp
c     write(*,*)'DOS',i,' met'
      i0=mxl*(i-1)+4
      read(12,'(//)')
   40 read(12,fmtx)
     &    (dmy,(dat(l,k,is),l=i0,i0+mxl-1),k=1,mse)
      endif
      if(index(a,'total DOS') .gt. 0) then
c     write(*,*)'total DOS met'
      read(12,'(1x,f12.7,3x,f13.5)')(dat(1,k,is),dat(2,k,is),k=1,mse-1)
      dat(1,mse,is)=dat(1,mse-1,is)
      dat(2,mse,is)=dat(2,mse-1,is)
      go to 20
      endif
   30 read(12,'(a)',end=10)a
   20 read(12,'(a)',end=10)a
      if(ismx .eq. 1) then
      do 60 i=1,mxl*icmp+3
      do 60 k=1,mse
   60 dat(i,k,2)=dat(i,k,1)
      endif
      close(12)
c     --- fix axis design
      xmi=big
      xmx=-big
      ymi=big
      ymx=-big
      do 70 is=1,2
      do 70 k=1,mse
      xmi=min(xmi,dat(1,k,is))
      xmx=max(xmx,dat(1,k,is))
      ymi=min(ymi,dat(2,k,is))
   70 ymx=max(ymx,dat(2,k,is))
      ymi=0d0
c     call ezaxis(xmi,xmx,nex,ntickx,1d0)
      call ezaxis(ymi,ymx,ney,nticky,0d0)
c     xmi=xmi*10d0**nex
c     xmx=xmx*10d0**nex
c     xincr=(xmx-xmi)/dble(ntickx-1)
      ymi=ymi*10d0**ney
      ymx=ymx*10d0**ney
      yincr=(ymx-ymi)/dble(nticky-1)
c     --- construct a gnuplot program
      open(13,file=file(1:len)//'.plt',status='unknown')
c     write(13,'(a)')'#!/usr/bin/gnuplot -persist'
      write(13,'(a)')'#!/usr/bin/gnuplot'
      write(13,'(2a)')'# set terminal postscript landscape noenhanced ',
     &               'monochrome dashed defaultplex "Helvetica" 14'
      write(13,'(a)')'# set output "dos.eps"'
      write(13,'(a)')'set term x11'
      write(13,'(a)')'set yzeroaxis lt -1 lw 0.1'
c     write(13,'(a)')'set yzeroaxis lt -1 lw 0.1'
      write(13,'(a)')'set border 15 lw 0.1'
      write(13,'(a)')'set mxtics 2'
      write(13,'(a)')'set mytics 2'
c     write(13,'(a,e15.7,a,e15.7)')'set xtics',xmi,',',xincr
      write(13,'(a,e15.7,a,e15.7)')'set ytics',ymi,',',yincr
      write(13,'(a)')'set size 1,1'
      write(13,'(a)')'set origin 0,0'
      write(13,'(a)')'set multiplot'
      write(13,'(a)')'set size 0.9,0.4'
c     write(13,'(a)')'set border 14 lw 0.1'
c     write(13,'(a)')'set lmargin at screen 0.1'
c     write(13,'(a)')'set rmargin at screen 0.9'
      write(13,'(a)')'set origin 0.1,0.5'
      write(13,'(a)')'set lmargin 0'
      write(13,'(a)')'set rmargin 8'
      write(13,'(a)')'set tmargin 0'
      write(13,'(a)')'set bmargin 0'
      write(13,'(3a)')'set title "',file(1:len),'" font "arial,14"'
      write(13,'(a)')'set format x ""'
c     write(13,'(a,e15.7,a,e15.7,a)')'set xrange [',xmi,':',xmx,']' 
      write(13,'(a,e15.7,a,e15.7,a)')'set xrange [*:*]' 
      write(13,'(a,e15.7,a,e15.7,a)')'set yrange [',ymi,':',ymx,']' 
      write(13,'(a)')'plot "-" using 1:2 title "" with lines linetype 1'
      write(13,fmt)((dat(l,k,1),l=1,mxl*icmp+3),k=1,mse)
      write(13,'(a)')'end'
      write(13,'(a)')'set border 11 lw 0.1'
      write(13,'(a)')'set origin 0.1,0.1'
      write(13,'(a)')'unset title'
      write(13,'(2a)')
     & 'set xlabel "Energy relative to Fermi energy (Ry)"',
     & ' font "arial,14"'
      write(13,'(2a)')
     & 'set label "DOS (1/Ry/unit cell)" at screen 0.04, screen 0.5',
     & ' center rotate font "arial,14"'
      write(13,'(a,e15.7,a,e15.7,a)')'set yrange [',ymx,':',ymi,']' 
      write(13,'(a)')'set format x'
      write(13,'(a)')'plot "-" using 1:2 title "" with lines linetype 1'
      write(13,fmt)((dat(l,k,2),l=1,mxl*icmp+3),k=1,mse)
      write(13,'(a)')'end'
      write(13,'(a)')'unset multiplot'
      write(13,'(a)')'pause -1'
c     write(13,'(a)')'#    EOF'
      close(13)
      call system('gnuplot '//file(1:len)//'.plt')
      deallocate(dat,stat=ierr)
      stop
   10 write(*,'(a/a)')'Illegal file specified','Usage: gpd file'
      close(12)
      end
      subroutine ezaxis(xmin,xmax,nex,ntick,adjust)
c-----------------------------------------------------------------------
c     I fix a scale of the axis.
c     In output, xmin is maximum integer for which  xmin*10**nex 
c     does not exeed original input xmin.
c     Similary, xmax is maximum integer for which  xmin*10**nex 
c     is not smaller than original xmax.
c     ntick gives a number of ticks including both end points.
c     Coded by H.Akai, April 28, 1991
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension nw(12),nh(12)
      data nw/ 4, 5, 6, 7, 8,10,12,14,16,20,25,30/
     &    ,nh/ 4, 5, 6, 7, 4, 5, 6, 7, 8, 4, 5, 6/
     &    ,small/1d-3/
c---- I fix scale of the axis.
c---- in some cases containing x=0 axis will give better presntation.
      if(xmax .lt. xmin) then
      write(6,1000)
 1000 format('   ***err in ezaxis...xmin and xmax not consistent')
      stop
      endif
      if(xmax*xmin .gt. 0d0 .and. 
     &    (xmax-xmin)/max(abs(xmax),abs(xmin)) .gt. adjust) then
      xmin=min(xmin,0d0)
      xmax=max(xmax,0d0)
      endif
      wx=xmax-xmin
c---- wx indicates x-data spreading.
      nex=int(log10(wx/2.9d0)+100d0+small)-100
      sx=10d0**nex
      xmin0=xmin/sx
      xmax0=xmax/sx
c---- now xmin and xmax is normalized by 10**nex so as to 
c---- satisfy 2.9<xmax-xmin<29.
      xmin=dble(int(xmin0+small))
      xmax=dble(int(xmax0-small))
c---- bug easily comes in the following statement. never
c---- use xmin instead of xmin0 in the logical-if, otherwise
c---- the procedure will fail when -1<xmin<0.
c---- xmin is maximum integer for which xmin*10**nex does not exeed
c---- original xmin.
      if(xmin0 .lt. 0d0) xmin=xmin-1d0
c---- xmax is maximum integer for which xmax*10**nex is not smaller
c---- than original xmax.
      if(xmax0 .gt. 0d0) xmax=xmax+1d0
      wx=xmax-xmin
c---- I then choose appropriate design for tick construction.
      iwx=int(wx+small)
      do 30 i=1,12
      ii=i
      if(nw(i) .ge. iwx) go to 40
   30 continue
c---- i adjust positioning such that the center of the data spreading
c---- is not much distant from that of the axis.
c---- however, if one of xmin or xmax is located near zero, xmin=0 or 
c---- xmax=0 is the arrangement of the highest priority.
   40 ws=dble(nw(ii))
      wm=ws-wx
      if(xmin .gt. -small .and. xmin .lt. wm+small) then
        wm=xmin
c---- xmin=0 is realized in this case.
      else if(xmax .lt. small .and. xmax .gt. -wm-small) then
        wm=xmax+wm
c---- xmax=0 is realized in this case.
      else
        if(xmax-xmax0 .gt. xmin0-xmin) wm=wm+1d0
        wm=dble(int(wm/2d0+small))
c---- otherwise, ticks are arranged so as to spread more or less 
c---- on centered.
      endif
      xmin=xmin-wm
      xmax=xmin+ws
      ntick=nh(ii)+1
c---- the first and the last ticks of axis are fixed.
c---- now I can start plotting. 
      return
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
