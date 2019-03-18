      subroutine mdexp2(da,e,l,aa,bb,s,n)
c-----------------------------------------------------------------------
c     Gives n-th stage of refinement of an extended midpoint rule.
c     See Numerical Recipes 2nd edition p135.
c     Modified by H. Akai, 13 April 1996, osaka
c  (C) Copr. 1986-92 Numerical Recipes Software !+!).
c     Modified to adapt OpenMP version by H. Akai, 6 Jan 2016, Tokyo.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dda=da**2
      dmy=bb
      b=exp(-aa)
      a=0d0
      if (n .eq. 1) then
      x=0.5d0*(a+b)
      y=-log(x)
      yy=y**2
      arg=-yy*dda+e/yy
      s=(b-a)*yy**l*exp(arg)/x
      else
      it=3**(n-2)
      tnm=it
      del=(b-a)/(3d0*tnm)
      ddel=del+del
      x=a+0.5d0*del
      sum=0d0
      do 10 j=1,it
      y=-log(x)
      yy=y**2
      arg=-yy*dda+e/yy
      if(arg .gt. -60d0) sum=sum+yy**l*exp(arg)/x
      x=x+ddel
      y=-log(x)
      yy=y**2
      arg=-yy*dda+e/yy
      if(arg .gt. -60d0) sum=sum+yy**l*exp(arg)/x
      x=x+del
   10 continue
      s=(s+(b-a)*sum/tnm)/3d0
      endif
      return
      end
