      subroutine rotatm(a,b,u,n)
c-----------------------------------------------------------------------
c     rotate atom following the operation u.
c     Ceded by H.Akai, 24 Dec. 96, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(3,n),b(3,n),u(3,3)
c     write(*,'((1x,3f10.5))')u
c     write(*,'((1x,9f8.5))')((a(i,j),i=1,3),j=1,n)
      do 10 i=1,n
      b(2,i)=u(1,1)*a(2,i)+u(1,2)*a(3,i)+u(1,3)*a(1,i)
      b(3,i)=u(2,1)*a(2,i)+u(2,2)*a(3,i)+u(2,3)*a(1,i)
   10 b(1,i)=u(3,1)*a(2,i)+u(3,2)*a(3,i)+u(3,3)*a(1,i)
      end
