      subroutine diffn(f,fp,  xr, dr, n )
c-----------------------------------------------------------------------
c      differentiate function f
c      Hubert Ebert, Muenchen
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension f(n),fp(n),xr(n),dr(n)
      nm2=n-2
c
c     forward difference at the beginning of the table
c
      fp(1) = ((2d0*f(4)+18d0*f(2))-(9d0*f(3)+11d0*f(1)))/6d0
      fp(2) = ((2d0*f(5)+18d0*f(3))-(9d0*f(4)+11d0*f(2)))/6d0
c
c     central difference at the interior of the table
c
      do 10 i=3,nm2
      fp(i)=((f(i-2)+8d0*f(i+1))-(8d0*f(i-1)+f(i+2)))/12d0
   10 continue
c
c     backward difference at the end of the table
c
      fp(n)  =((11d0*f(n  )+9d0*f(n-2))-(18d0*f(n-1)+2d0*f(n-3)))/6d0
      fp(n-1)=((11d0*f(n-1)+9d0*f(n-3))-(18d0*f(n-2)+2d0*f(n-4)))/6d0
      do 20 i=1,n
   20 fp(i) = fp(i) / dr(i)
      end
