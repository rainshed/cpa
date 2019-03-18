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
