      function isrtyp(reltyp)
c-----------------------------------------------------------------------
c     Given a reltyp name, this program returns the corresponding code.
c     coded by H.Akai, April 1992, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(ityp=4)
      character*5 relval(ityp),reltyp
      data relval/'nrl','sra','nrlls','srals'/
      do 10 i=1,ityp
      isrtyp=i-1
      if(reltyp .eq. relval(i)) return
   10 continue
      call errtrp(2,'isrtyp','illegal reltyp: nrl assumed')
      isrtyp=0
      end
