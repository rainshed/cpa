      subroutine rlubks(a,n,indx,b)
c-----------------------------------------------------------------------
c     real*8 version of lu backward substitution
c  (C) Copr. 1986-92 Numerical Recipes Software !+!).
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer indx(n)
      real*8 a(n,n),b(n)
      ii=0
      do 20 i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii .ne. 0) then
      do 10 j=ii,i-1
      sum=sum-a(i,j)*b(j)
   10 continue
      else if (abs(sum) .gt. 0d0) then
      ii=i
      endif
      b(i)=sum
   20 continue
      do 40 i=n,1,-1
      sum=b(i)
      do 30 j=i+1,n
      sum=sum-a(i,j)*b(j)
   30 continue
      b(i)=sum/a(i,i)
   40 continue
      return
      end
