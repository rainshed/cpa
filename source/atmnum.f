      function atmnum(atom)
c----------------------------------------------------------------------
c     Inverse function of asymbl
c     coded by H.Akai, 1985, Juelich
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character atom*2,asymbl*2
      atmnum=0d0
      do 10 i=0,104
      z=dble(i)
      if(atom .eq. asymbl(z)) then
      atmnum=z
      return
      endif
   10 continue
      return
      end
