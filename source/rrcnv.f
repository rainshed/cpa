      subroutine rrcnv(ro,msr,meshr)
c-----------------------------------------------------------------------
c     r(up,down) --> r(charge,spin) conversion
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ro(msr,2)
      do 10 k=1,meshr
      s=ro(k,1)+ro(k,2)
      ro(k,2)=ro(k,1)-ro(k,2)
   10 ro(k,1)=s
      return
      end
