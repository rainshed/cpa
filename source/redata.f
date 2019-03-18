      function redata(buff)
c-----------------------------------------------------------------------
c     Given character data 'buff', this program convert the data into
c     numerical data.
c     The character data could have the following forms:
c
c     1) Fixed and floating data such as 1, 1d0, 1.0, etc.
c     2) Fractional data such as 1/3, 1d0/3d0, 2/3d0, etc.
c
c     Coded by H. Akai, 22 July 2007, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character buff*(*)
      n=len(buff)
      isl=0
      do 10 i=1,n
      if(buff(i:i) .eq. '/') then
      buff(i:i)=' '
      isl=isl+1
      endif
   10 continue
      if(isl .eq. 0) then
      read(buff,*)redata
      else if(isl .eq. 1) then
      read(buff,*)a,b
      redata=a/b
      else
      call errtrp(1,'redata','illegal format found in buffer')
      endif
      end
