      subroutine laguer(a,m,x)
c-----------------------------------------------------------------------
c     Laguere algorithm of solving all zero's of a polinomial.
c     See Numerical Recipes.
c     modified and addapted to KKR packege by H.Akai, 1993
c     Revised by H. Akai, Tokyo, Sep. 2013
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (epss=3d-16,mr=8,mt=10,maxitr=mt*mr)
      complex*16 x,xini,dx,x1,b,d,f,g,h,sq,gp,gm,g2
      real*8 frac(mr),a(m+1)
      data frac /.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1d0/
     &    ,small/1d-20/, inimx/100/, xincr/1d-1/
      dm=dble(m)
c     write(*,'(1x,a/(1x,1p,6e13.6))')
c    &    'polinomial coefficients',(a(k),k=1,m)
      xini=x
      do 30 ini=1,inimx
c     --- in the case that the iteration does not converge,
c         start from a different trial x.
      x=xini+xincr*dble(ini-1)
      do 20 itr=1,maxitr
c     write(*,'(1x,a,i3,a,1p,2e13.6)') 'itr=',itr,'  guess=',x
      b=dcmplx(a(m+1),0d0)
      err=abs(b)
      d=dcmplx(0d0,0d0)
      f=dcmplx(0d0,0d0)
      abx=abs(x)
      do 10 j=m,1,-1
      f=x*f+d
      d=x*d+b
      b=x*b+a(j)
   10 err=abs(b)+abx*err
      err=epss*err
      if(abs(b).le.err) then
      return
      else
      g=d/b
      g2=g*g
      h=g2-2d0*f/b
      sq=sqrt((dm-1d0)*(dm*h-g2))
      gp=g+sq
      gm=g-sq
      abp=abs(gp)
      abm=abs(gm)
      if(abp.lt.abm) gp=gm
      if (max(abp,abm).gt.0d0) then
      dx=dm/gp
      else
      dx=exp(dcmplx(dlog(1d0+abx),dble(itr)))
      endif
      endif
      x1=x-dx
      if(abs(x-x1) .lt. small)return
      if (mod(itr,mt).ne.0) then
      x=x1
      else
      x=x-dx*dble(frac(itr/mt))
      endif
   20 continue
   30 continue
      call errtrp(1,'laguer','no convergence')
      end
