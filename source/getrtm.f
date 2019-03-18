      subroutine getrtm(u,alpha,beta,gamma,mxl)
c-----------------------------------------------------------------------
c     gets transformation matrix of cubic operation for real harmonics.
c     coded by H.Akai, 14 Dec. 1996, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mmxl=5)
      real*8 u((2*mxl-1)**2,mxl),v1(3,2*mmxl-1),v2(3),p((2*mmxl-1)**2)
     &      ,a(mmxl**2,2*mmxl-1),b(mmxl**2,2*mmxl-1)
      logical init
      save init,v1,mxmg
      data init/.true./
      if(mxl .gt. mmxl) call errtrp(1,'getrtm','mxl too large')
      if(init) then
      init=.false.
      call getnum(1996,v1,3*(2*mmxl-1))
      endif
      mxmg=2*mxl-1
      do 10 i=1,mxmg
      call vrotat(v2,v1(1,i),alpha,beta,gamma)
      call realh(v1(1,i),a(1,i),mxl-1)
   10 call realh(v2,b(1,i),mxl-1)
      do 20 l=1,mxl
      n=2*l-1
      i=0
      do 30 k1=(l-1)**2+1,l**2
      do 30 k2=1,n
      i=i+1
      p(i)=a(k1,k2)
   30 u(i,l)=b(k1,k2)
      call drvlud(p,u(1,l),n)
   20 call transp(u(1,l),n)
      end
