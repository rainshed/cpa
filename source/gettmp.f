      subroutine gettmp(lunit)
c----------------------------------------------------------------------
c     This program assign temporary file to unit=luit.
c     coded by H.Akai, Dec. 1992, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      close(lunit)
      open(lunit,form='unformatted',status='scratch')
      return
      end
