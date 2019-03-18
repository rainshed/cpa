      function ifkey(key,keys)
c-----------------------------------------------------------------------
c     Given key and keys this program returns if string key
c     is involved in string keys. If it is the case, key
c     is removed from keys.
c
c     *** in this version the key is not removed.
c     Coded by H. Akai, 23 Nov. 2007
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character keys*(*),key*(*),dmys*100,dmy*100
      logical ifkey
      ifkey=.false.
      dmys=keys
      dmy=key
      m=lftfil(dmys)
      n=lftfil(dmy)
      i=index(dmys(1:m),dmy(1:n))
      if(i .gt. 0) then
      ifkey=.true. 
c     keys(i:i+n-1)=' '
c     m=lftfil(keys)
      endif
      end
