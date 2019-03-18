      function cg(l1,l2,l3,m1,m2,m3)
c-----------------------------------------------------------------------
c     Given l1,l2,l3,m1,m2,m3, this function returns the gount number.
c     coded by M.Akai 1979, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      s(n1)=sign(1d0 , dble(n1)+1d-1)
      t(n1)=s(n1)*dble( (-1)**n1 )
      u(n1,n2)=-dble( (-1)**(n1+n2) )/sqrt(2d0)
      pi=4d0*atan(1d0)
      cg=0d0
      if(l1.ge.iabs(m1).and.l2.ge.iabs(m2).and.l3.ge.iabs(m3)) go to 10
      return
c
   10 if(l3.lt.iabs(l1-l2).or.l3.gt.l1+l2 .or.s(m1)*s(m2)*s(m3).lt.0d0)
     &  return
c
      n1=m1
      n2=m2
      n3=m3
      if(m1+m2 .eq. m3) go to 20
      if(m1+m2 .eq. -m3) go to 30
      if(m1-m2 .eq. m3) go to 40
      if(m1-m2 .eq. -m3) go to 50
      return
c
   20 cg=1d0
      go to 60
c*****
   30 cg=t(m3)
      n3=-n3
      go to 60
c*****
   40 cg=t(m2)
      n2=-n2
      go to 60
c*****
   50 cg=t(m1)
      n2=-n2
      n3=-n3
   60 if(m1 .ge. 0 .and. m2 .ge. 0) go to 70
      if(m1 .ge. 0) go to 80
      if(m2 .ge. 0) go to 90
      cg=cg*u(m1,m2)
      go to 100
c*****
   70 cg=cg/sqrt(2d0)
      go to 100
c*****
   80 cg=cg*u(m2,m3)
      go to 100
c*****
   90 cg=cg*u(m1,m3)
  100 if(m1 .eq. 0) cg=cg*sqrt(2d0)
      if(m2 .eq. 0) cg=cg*sqrt(2d0)
      if(m3 .eq. 0) cg=cg*sqrt(2d0)
      if(m1 .eq. 0 .and. m2 .eq. 0 .and. m3 .eq. 0) cg=cg/2d0
      cg=cg*sqrt( dble((2*l1+1)*(2*l2+1))/dble(2*l3+1)/4d0/pi)
     &    *s(m3)*cgc(l1,l2,l3,n1,n2,n3)*cgc(l1,l2,l3,0,0,0)
      return
      end
