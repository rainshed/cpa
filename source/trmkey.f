      function trmkey(key,keys)
c-----------------------------------------------------------------------
c     Given key and keys this program chek if string key
c     is involved in string keys. If it is the case, trmkey
c     return the value of keys after removing key from keys.
c     Coded by H. Akai, Tokyo, 26 Jan. 2017.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character keys*(*),key*(*),trmkey*(*),dmys*100,dmy*100
      logical ifkey
      trmkey=keys
      dmy=key
      m=lftfil(trmkey)
      n=lftfil(dmy)
      i=index(trmkey(1:m),dmy(1:n))
      if(i .gt. 0) then
      trmkey(i:i+n-1)=' '
      m=lftfil(trmkey)
      endif
      end
