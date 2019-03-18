      subroutine radnrl(e,j,r,kmatch,g1,g2,node,z,a,b,xr,msr,ns1)
c-----------------------------------------------------------------------
c     This program calculaets the non-relativistic radial wave function
c     for a given potential and an energy by the Noumerov method.
c      e      ... energy
c      j      ... =l+1
c      r      ... =y*x/sqrt(x'), where y is the radial wave function,
c                 x(s) is any  mesh point specified by s, x' indicates
c                 its derivative with respect to s.  r is not normalized
c      kmatch ... Indicating the mesh point where the outward and the
c                 inward integrations meet. if kmatch < or = 0, a
c                 suitable kmatch is generated in the program. otherwise
c                 a given one is used.
c      g1, g2 ... Logarithmic derivatice of r(outward) and r(inward),
c                 with respect to s, s indicating the mesh points,
c                 point.
c      node   ... The number of the nodes of the wave function.
c      z      ... Specifying the potential v(x)=-z(x)/x
c      a,b,xr ... Radial mesh must have a form xr(s)=a*(exp(b*s)-1)
c      msr    ... s=1,2,3, ..., msr
c      ns1    ... Array size of r, z and xr.
c
c     With the transformation x ---> s and y ---> r indicated above,
c     the Schroedinger eauation becomes  r''(s)=b**2*gm(s)*r(s), where
c     gm(s) is the function given by the statement function in the
c     program. The procedure corresponding to the numerov method
c     is also given in the form of a statement function.
c
c     This program is not optimized for array processor.
c
c             coded by H.Akai  on Jan. 8, 1986  (Osaka)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension r(ns1),z(ns1),xr(ns1),str(3)
      logical chk,ok
      data bigexp/8d1/
c
c     ----- statement function -----
      gm(arg1,arg2)=2.5d-1+(arg1+a)**2*((-arg2+al/arg1)/arg1-e)
      rnxt(arg1,arg2)=((bb24+10d0*gm2)*arg2-(bb12-gm1)*arg1)/(bb12-gm3)
c     ------------------------------
c
      arg1=0d0
      arg2=0d0
      l=j-1
      al=dble(l*j)
      bb12=12d0/b**2
      bb24=2d0*bb12
      zn=5d-1*z(1)
c     --- the following statement to cheat compilers ---
      zn=zn+arg1+arg2
c
c      --- emx is a number for which exp(-sqrt(emx)) is very small ---
      emx=bigexp**2*zn**(-.66666666666d0)
c
c     --- power expansion near x=0 ---
c        y=x**l*(1-z*r/(l+1)+ ... )
c     the following expansion corresponds to y ---> r transformation
c     of the above expression.
c
      r(1)=0d0
      r(2)=xr(2)**j*(1d0-(zn/dble(j)+5d-1/a)*xr(2))
      gm1=gm(xr(2),z(2))
      gm2=gm(xr(3),z(3))
      r(3)=(bb24+10d0*gm1)*r(2)
      if(j .eq. 1) r(3)=r(3)-2d0*zn*a**2
      r(3)=r(3)/(bb12-gm2)
c
c     --- outward integration by the noumerov method ---
      ok=.false.
      is1=1
      node=0
      do 10 k=4,msr-2
      kk=k
      gm3=gm(xr(k),z(k))
      r(k)=rnxt(r(k-2),r(k-1))
      gm1=gm2
      gm2=gm3
      if(kmatch .ge. 4) go to 20
      chk=gm3 .gt. 0d0
      if(.not. chk) ok=.true.
      if(chk .and. ok) go to 30
   20 if(k .eq. kmatch) go to 30
      is2=sign(1.1d0,r(k))
      if(is1 .ne. is2) node=node+1
   10 is1=is2
   30 kmatch=kk
      do 40 k=kk+1,kk+2
      gm3=gm(xr(k),z(k))
      r(k)=rnxt(r(k-2),r(k-1))
      gm1=gm2
   40 gm2=gm3
      do 50 k=1,3
   50 str(k)=r(kk+k-3)
      g1=(r(kk-2)-8d0*(r(kk-1)-r(kk+1))-r(kk+2))/12d0/r(kk)
      do 60 k=1,msr-kmatch-1
      kst=k
      kk=msr-k+1
      r(kk)=0d0
      if(al-(e*xr(kk)+z(kk))*xr(kk) .lt. emx) go to 70
   60 continue
   70 r(kk-1)=1d-35
      klst=kk-1
c
c     --- inward integration by the noumerov method ---
      gm1=gm(xr(kk),z(kk))
      gm2=gm(xr(kk-1),z(kk-1))
      do 80 k=kst+2,msr-kmatch+3
      kk=msr-k+1
      gm3=gm(xr(kk),z(kk))
      r(kk)=rnxt(r(kk+2),r(kk+1))
      gm1=gm2
   80 gm2=gm3
      kk=kmatch
      g2=(r(kk-2)-8d0*(r(kk-1)-r(kk+1))-r(kk+2))/12d0/r(kk)
      red=str(3)/r(kk)
      do 90 k=1,3
   90 r(kk+k-3)=str(k)
      do 100 k=kk+1,klst
  100 r(k)=r(k)*red
      return
      end
