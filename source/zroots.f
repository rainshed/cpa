      subroutine zroots(a,m,roots,ierr)
c-----------------------------------------------------------------------
c     Given m and the real coefficients of a polinomial
c     a(1)+a(2)*x+a(3)*x**2+...+a(m+1)*x**m, this routine returns
c     all the roots of the polinomials. In output, nr is the
c     number of roots. It neccessarily is not equal to m since
c     a(m+1), etc. could be zero.
c     See Numerical Recipes (2nd ed.) pp365-368.
c     coded by H.Akai, 1993
c     Modified by H. Akai, Tokyo, Jan. 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (small=1d-7, zero=1d-10)
      real*8 a(m+1)
      real*8,allocatable::ad(:)
      complex*16 roots(m),x,z
      allocate(ad(m+1))
c     write(*,'(1x,a/(1x,1p,6e13.6))')
c    &    'polinomial coefficients',(a(k),k=1,m+1)
      do 10 i=1,m+1
   10 ad(i)=a(i)
      j=m
      do 20 i=1,m
      x=dcmplx(0d0,0d0)
      call laguer(ad,j,x)
      if(abs(dimag(x)) .le. zero*abs(dble(x))) then
      rx=dble(x)
      roots(j)=dcmplx(rx,0d0)
      b=ad(j+1)
      do 30 jj=j,1,-1
      c=ad(jj)
      ad(jj)=b
   30 b=rx*b+c
      j=j-1
      else
      roots(j)=x
      roots(j-1)=conjg(x)
      p=-2d0*dble(x)
      q=abs(x)**2
      b=ad(j+1)
      c=ad(j)
      do 70 jj=j,2,-1
      d=ad(jj-1)
      ad(jj-1)=b
      b=c-p*ad(jj-1)
   70 c=d-q*ad(jj-1)
      j=j-2
      endif
      if(j .le. 0) go to 80
   20 continue
   80 continue
c     write(*,*)'check zero for set 1'
c     do 24 j=1,m
c     z=a(m+1)
c     do 22 i=m,1,-1
c  22 z=z*roots(j)+a(i)
c  24 write(*,'(1x,i3,2e15.7)')j,z
      ierr=0
      do 40 j=1,m
      x=roots(j)
      call laguer(a,m,roots(j))
c     write(*,'(a,i2,4e15.7)')' j,old,new=',j,x,roots(j)
   40 if(abs(x-roots(j))/abs(x) .gt. small) ierr=1
c     write(*,*)'check zero for set 2'
c     do 64 j=1,m
c     z=a(m+1)
c     do 62 i=m,1,-1
c  62 z=z*roots(j)+a(i)
c  64 write(*,'(1x,i3,2e15.7)')j,z
      do 60 j=2,m
      x=roots(j)
      do 50 i=j-1,1,-1
      if(abs(roots(i)) .le. abs(x)) go to 60
   50 roots(i+1)=roots(i)
      i=0
   60 roots(i+1)=x
c     write(*,*)(roots(i),i=1,m)
      deallocate(ad)
      end
