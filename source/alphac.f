      function alphac(a)
c----------------------------------------------------------------------
c     Given a character, this program tells whether it contains only
c     alpha cahracters (i.e. A-Z, a-z, and space) or others.
c     Coded by H. Akai, 24 Oct. 98, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical alphac
      character a*(*)
      n=len(a)
      alphac=.false.
      do 10 i=1,n
      if((a(i:i) .lt. 'a' .or. a(i:i) .gt. 'z') .and.
     &    (a(i:i) .lt. 'A' .or. a(i:i) .gt. 'Z') .and.
     &    (a(i:i) .ne. ' ')) return
   10 continue
      alphac=.true.
      end
