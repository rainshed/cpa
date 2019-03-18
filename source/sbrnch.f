      subroutine sbrnch(a,m,k)
c-----------------------------------------------------------------------
c     Select a proper branch so as to keep continuity.
c     coded by H.Akai, 1986, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(m,k)
      pi=4d0*atan(1d0)
      do 10 i=1,m
      jump=int((dimag(a(i,1))+pi)/pi/2d0+1d5+1d-5)
      jump=jump-100000
      a(i,1)=a(i,1)-dble(jump)*2d0*pi*(0d0,1d0)
      do 10 j=2,k
      jump=int((dimag(a(i,j)-a(i,j-1))+pi)/pi/2d0+1d5+1d-5)
      jump=jump-100000
c     if(jump .ne. 0) write(*,'(1x,2i4,2f12.5)')
c    & i,j,dimag(a(i,j))/pi,dimag(a(i,j-1))/pi
   10 a(i,j)=a(i,j)-dble(jump)*2d0*pi*(0d0,1d0)
c  10 a(i,j)=a(i,j)
      return
      end
