      subroutine gpsort(gpt,ngpt)
c-----------------------------------------------------------------------
c     This program sort the g vector in ascending order.
c     coded by H.Akai, April 1992, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 gpt(3,ngpt),b(3)
      real*8,allocatable::d(:)
      allocate(d(ngpt))
      do 10 i=1,ngpt
   10 d(i)=gpt(1,i)**2+gpt(2,i)**2+gpt(3,i)**2
c     --- sort by length of g vector
c     --- simple sorting is most efficient for such a small problem.
      do 20 j=2,ngpt
      a=d(j)
      do 30 l=1,3
   30 b(l)=gpt(l,j)
      do 40 i=j-1,1,-1
      if(d(i) .lt. a) go to 50
      d(i+1)=d(i)
      do 60 l=1,3
   60 gpt(l,i+1)=gpt(l,i)
   40 continue
      i=0
   50 d(i+1)=a
      do 70 l=1,3
   70 gpt(l,i+1)=b(l)
   20 continue
      deallocate(d)
      end
