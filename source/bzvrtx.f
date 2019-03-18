      subroutine bzvrtx(g,vx,nv)
c-----------------------------------------------------------------------
c     Given a g-vector, this program constructs the first Brilloiun-zone
c     and returns the positions ((v(i,n),i=1,3),n=1,nv) of all the
c     vertices of BZ polihedron. nv is the number of the vortices.
c     Coded by H. Akai, 3 Dec. 2014, Tokyo.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(npx=27,nvx=24)
c     --- npx=27 is enough. nv never exceeds 24. nvx=24 or larger
c         should be specified in the calling routines.
      real*8 p(3,npx),d(npx),a(3,3),b(3),vx(3,nvx),g(3,3)
      integer indx(3)
      data tiny/1d-12/
      np=0
      do 10 l=-1,1
      do 10 m=-1,1
      do 10 n=-1,1
      if(l**2+m**2+n**2 .ne. 0) then
      np=np+1
      do 20 i=1,3
   20 p(i,np)=dble(l)*g(i,1)+dble(m)*g(i,2)+dble(n)*g(i,3)
      d(np)=5d-1*(p(1,np)**2+p(2,np)**2+p(3,np)**2)
      endif
   10 continue
      nv=0
      do 40 l=1,np-2
      do 40 m=l+1,np-1
      do 40 n=m+1,np
c     --- three planes are chosen.
      do 50 i=1,3
      a(1,i)=p(i,l)
      a(2,i)=p(i,m)
   50 a(3,i)=p(i,n)
c     --- position of the vertex formed by the three planes are obtained.
      b(1)=d(l)
      b(2)=d(m)
      b(3)=d(n)
      call rludcm(a,3,indx,dd)
      if(abs(dd*a(1,1)*a(2,2)*a(3,3)) .lt. tiny) go to 40
      call rlubks(a,3,indx,b)
c     --- check if the vertex and (0,0,0) is on the same side of all
c         other planes.
      do 60 j=1,np
      if((j-l)*(j-m)*(j-n) .eq. 0) go to 60
      if(d(j)*(d(j)-p(1,j)*b(1)-p(2,j)*b(2)
     &   -p(3,j)*b(3)) .lt. -tiny) go to 40
   60 continue
      do 70 j=1,nv
      if((b(1)-vx(1,j))**2+(b(2)-vx(2,j))**2+(b(3)-vx(3,j))**2
     &   .lt. tiny) go to 40
   70 continue
      nv=nv+1
      if(nv .gt. nvx) call errtrp(1,'bzvrtx','nv too large')
      do 80 i=1,3
   80 vx(i,nv)=b(i)
   40 continue
c    --- now all the vertices are found.
c     write(*,'(a,i3)')'number of vrtices =',nv
c     write(*,'(i3,3e13.5)')(n,(vx(i,n),i=1,3),n=1,nv)
      end
