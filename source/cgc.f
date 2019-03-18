      function cgc(l1,l2,l3,m1,m2,m3)
c-----------------------------------------------------------------------
c     Clebsh-Gordan coefficient covering l1+l2+l3 <= lsmx.
c     For bigger lsmx, nel also must be increased,
c     e.g., nel=22 can cover factorial of 81.
c     coded by M.Akai and H.Akai, 1979, Osaka
c     revised by H.Akai 29 Dec. 1990
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (nel=22, lsmx=80, ntb=lsmx+2)
      integer n(nel,ntb),l(nel),k(nel),ipr(nel)
      logical init
      save init,ipr,n
      data init/.true./
      if(l1+l2+l3 .gt. lsmx)
     &  call errtrp(1,'cgc','l1+l2+l3 is too large')
c     --------------------------------------------------------
c     this part runs once for all, constructing the table.
c     --------------------------------------------------------
      if(init)then
      init=.false.
c
c     get prime numbers
      jel=1
      ipr(1)=2
      do 30 i=3,lsmx+1
      imx=sqrt(dble(i)+5d-1)
      do 10 m=1,jel
      if(ipr(m) .gt. imx) go to 20
      if(mod(i,ipr(m)) .eq. 0) go to 30
   10 continue
   20 jel=jel+1
      if(jel .gt. nel) call errtrp(1,'cgc','nel too small')
      ipr(jel)=i
   30 continue
c
c     foctorizing factorials
      do 40 i=1,jel
   40 n(i,1)=0
      do 60 j=2,ntb
      jj=j-1
      do 60 i=1,jel
      n(i,j)=n(i,j-1)
      do 50 m=1,10000
      if(mod(jj,ipr(i)) .ne. 0) go to 60
      jj=jj/ipr(i)
   50 n(i,j)=n(i,j)+1
   60 continue
      endif
c     --------------------------------------------------------
c
      cgc=0d0
      if(l1 .lt. 0 .or. l2 .lt. 0 .or. l3 .lt. 0) return
      if(l3 .lt. abs(l1-l2) .or. l3 .gt. l1+l2) return
      if(abs(m1) .gt. l1 .or. abs(m2) .gt. l2
     &   .or. abs(m3) .gt. l3) return
      if(m1+m2-m3 .ne. 0) return
      if(mod(l1+l2+l3,2) .ne. 0) return
      jel=nel
      do 70 i=nel,2,-1
   70 if(ipr(i) .gt. l1+l2+l3+1) jel=i-1
      j=2*l3+2
      do 80 i=1,jel
   80 l(i)=n(i,j)-n(i,j-1)-n(i,l1+l2+l3+2)
     &     +n(i,l1+l2-l3+1)+n(i,l3+l1-l2+1)+n(i,l2+l3-l1+1)
     &     +n(i,l1+m1+1)+n(i,l1-m1+1)+n(i,l2+m2+1)
     &     +n(i,l2-m2+1)+n(i,l3+m3+1)+n(i,l3-m3+1)
      i1=max0(-l3+l2-m1, -l3+l1+m2, 0)+1
      i2=min0(l1+l2-l3, l1-m1, l2+m2)+1
      do 90 i=1,jel
   90 k(i)=0
      nm=0
      do 110 mm=i1,i2
      m=mm-1
      ns=(-1)**m
      do 100 i=1,jel
      jd=n(i,m+1)+n(i,l1+l2-l3-m+1)
     &   +n(i,l1-m1-m+1)+n(i,l2+m2-m+1)
     &   +n(i,l3-l2+m1+m+1)+n(i,l3-l1-m2+m+1)-k(i)
      if(jd .lt. 0) then
      ns=ns*ipr(i)**(-jd)
      else
      k(i)=k(i)+jd
      nm=nm*ipr(i)**jd
      endif
  100 continue
  110 nm=nm+ns
      prfns=1d0
      prfnl=1d0
      prfnk=1d0
      do 120 i=1,jel
      ls=l(i)/2
      prfns=prfns*dble(ipr(i))**ls
      prfnl=prfnl*dble(ipr(i))**(l(i)-2*ls)
  120 prfnk=prfnk*dble(ipr(i))**k(i)
      cgc=prfns*sqrt(prfnl)*dble(nm)/prfnk
      return
      end
