      function lwcase(token)
c----------------------------------------------------------------------
c     Given a character string token, this function returns its lower
c     case conversion if it contains upper case character.
c     coded by H.Akai, Feb. 1993, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character token*(*),lwcase*80
      ishift=ichar('a')-ichar('A')
      lwcase=token
      do 10 i=1,80
      if(lwcase(i:i) .eq. ' ') return
   10 if(lwcase(i:i) .ge. 'A' .and. lwcase(i:i) .le. 'Z')
     &   lwcase(i:i)=char(ichar(lwcase(i:i))+ishift)
      return
      end
