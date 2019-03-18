      subroutine midexp(funk,aa,bb,s,n)
c-----------------------------------------------------------------------
c     Gives n-th stage of refinement of an extended midpoint rule.
c     See Numerical Recipes 2nd edition p135.
c     Modified by H. Akai, 13 April 1996, osaka
c  (C) Copr. 1986-92 Numerical Recipes Software !+!).
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      external funk
      func(x)=funk(-log(x))/x
      dmy=bb
      b=exp(-aa)
      a=0d0
      if (n .eq. 1) then
      s=(b-a)*func(0.5d0*(a+b))
      else
      it=3**(n-2)
      tnm=it
      del=(b-a)/(3d0*tnm)
      ddel=del+del
      x=a+0.5d0*del
      sum=0d0
      do 10 j=1,it
      sum=sum+func(x)
      x=x+ddel
      sum=sum+func(x)
      x=x+del
   10 continue
      s=(s+(b-a)*sum/tnm)/3d0
      endif
      return
      end
