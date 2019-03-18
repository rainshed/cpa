      subroutine gengpt(g,gpt,ngptmx,ngpt,gmx,cornr,nc)
c----------------------------------------------------------------------
c     Given three primitive translation vectors, this program generates
c     a set of lattice vectors, whose maximum length do not exceeds
c     'rmx'. The vectors are stored in rpt in ascending order of the
c     length.
c     coded by H.Akai April 1992, Osaka
c     tiny modification by H.Akai, 6 Dec. 1995, Osaka
c     major modification by H.Akai, 1 Jul. 2007, Osaka
c     Modified so that the location of g vector preserve the symmetry
c     of crystal.
c     Modified by H. Akai, 4 Dec. 2014, Tokyo
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mx=50)
      real*8 g(3,3),gpt(3,ngptmx),v(3,3),cornr(3,nc)
c     Ggenerate the points for which the distance from all
c     the vertices is within gmx.
      ggmx=gmx**2
      ngpt=0
      do 10 l=-mx,mx
      do 20 i=1,3
   20 v(i,1)=dble(l)*g(i,1)
      do 10 m=-mx,mx
      do 30 i=1,3
   30 v(i,2)=v(i,1)+dble(m)*g(i,2)
      do 10 n=-mx,mx
      do 40 i=1,3
   40 v(i,3)=v(i,2)+dble(n)*g(i,3)
      vv=1d10
      do 70 k=1,nc
   70 vv=min(vv,(v(1,3)-cornr(1,k))**2+(v(2,3)-cornr(2,k))**2
     &           +(v(3,3)-cornr(3,k))**2)
      if(vv .lt. ggmx) then
      ngpt=ngpt+1
      if(ngpt .gt. ngptmx) call errtrp(1,'gengpt','ngpt too large')
      do 50 i=1,3
   50 gpt(i,ngpt)=v(i,3)
      endif
   10 continue
c     --- sort them in ascending order of the length
      call gpsort(gpt,ngpt)
      return
      end
