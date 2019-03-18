      subroutine cgtabl(ncub,clks,lmax,last)
c-----------------------------------------------------------------------
c     Generate the table used for the sum
c           a(l1,l2)=sum_l3 c(l1,l2,l3)*d(l3),
c     where c(l1,l2,l3) is the Gaunt number.
c     ----------------------------------------
c      lmax  1   2   3   4    5     6     7
c     ----------------------------------------
c      last 37  243 964 2854 6901 14723 28462
c     ----------------------------------------
c     coded by M.Akai, 1979, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ncub(*),clks(*)
      ll=(lmax+1)**2
      lj=(2*lmax+1)**2
      last=0
      do 20 mr=1,ll
      mrr=mr
      call subscr(l1,m1,mrr)
      do 20 mc=1,ll
      mcc=mc
      call subscr(l2,m2,mcc)
      do 10 mj=1,lj
      mjj=mj
      call subscr(l3,m3,mjj)
      c=cg(l1,l2,l3,m1,m2,m3)
      if(abs(c) .lt. 1d-7) go to 10
      last=last+1
      ncub(last)=mj
      clks(last)=c
   10 continue
      last=last+1
      ncub(last)=0
   20 clks(last)=0d0
c     write(*,*)' last=',last
      end
