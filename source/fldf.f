      function fldf(ro)
c-----------------------------------------------------------------------
c     Non-spin-polarized local density functional potential.
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      data
c     ---- vbh parametrization ----
c    &     c,a/.0504d0, 30d0/
c     ---- hl ( same as mjw ) parametrization ----
     &     c,a/.0449d0, 21d0/
      fldf=0d0
      if(ro .lt. 1d-20) return
      rs=(3d0/ro)**.333333333333d0
      fldf=-1.221774d0/rs-c*log(1d0+a/rs)
c
c     ---- slater's x-alpha with alpha=1 (herman-skillman) ----
c     fldf=-1.221774d0/rs*1.5d0
      return
      end
