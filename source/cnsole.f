      subroutine cnsole(iunit,buffer,*)
c----------------------------------------------------------------------
c     This program reads in character string from console (unit=5).
c     Lines started with 'c', 'c' or '#' are regarded as comment
c     lines and ignored. 'return 1' is executed when end of file is
c     detected.
c     Blank cards are simply ignored.
c     Commnet line rule should be changed in future. Many users
c     apt to use inputs starting by C. For example
c     to specify the name of site such as 'Cu'. This obviously
c     is understood as a comment line and lead an input error.
c     I think admitting only '#' as comment line would be safer.
c     coded by H.Akai, April 1992, Osaka
c     Last modified by H.Akai, 5 Aug. 1999, Duisburg
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character buffer*(*)
      buffer=' '
   10 read(iunit,'(a)',end=20)buffer
      do 30 i=1,len(buffer)
   30 if(ichar(buffer(i:i)) .lt. 32 .or.
     &   ichar(buffer(i:i)) .ge. 128) buffer(i:i)=' '
      if(buffer(1:1) .eq. 'c' .or. buffer(1:1) .eq. 'C'
     &  .or. buffer(1:1) .eq. '#' .or. buffer .eq. ' ') go to 10
      return
   20 return 1
      end
