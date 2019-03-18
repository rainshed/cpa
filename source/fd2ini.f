      function fd2ini(da,e,l)
c-----------------------------------------------------------------------
c     integrand for integral appearing in the Ewalt sum
c     Coded by H.Akai, 12 April 1996, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(zeroln=-60d0)
      save dda,ee,ll
c     initialization
      dda=da**2
      ee=e
      ll=l
      fd2ini=0d0
      return
c-----------------------------------------------------------------------
      entry fd2(x)
c-----------------------------------------------------------------------
      fd2=0d0
      xx=x**2
      arg=-xx*dda+ee/xx
      if(arg .lt. zeroln) return
      fd2=xx**ll*exp(arg)
      end
