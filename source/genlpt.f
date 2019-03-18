      subroutine genlpt(r,rpt,nrptmx,nrpt,rmx)
c----------------------------------------------------------------------
c     Given three primitive translation vectors, this program generates
c     a set of lattice vectors, whose maximum length do not exceeds
c     'rmx'. The vectors are stored in rpt in ascending order of the
c     length.
c     coded by H.Akai April 1992, Osaka
c     tiny modification by H.Akai, 6 Dec. 1995, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mx=50)
      dimension r(3,3),rpt(3,nrptmx),v(3,3),p(3,27)
c     First, find vertices of parallelepiped formed by vectors r.
c     The parallelepiped is located both  positive or negative
c     side of the cartesian cordinate.
      k=0
      do 60 l=-1,1
      do 60 m=-1,1
      do 60 n=-1,1
      k=k+1
      do 60 i=1,3
   60 p(i,k)=dble(l)*r(i,1)+dble(m)*r(i,2)+dble(n)*r(i,3)
c     write(*,'(1x,3f12.5)')p
c   
c     Then generate the points for which the distance from all
c     the vertices is within rmx.
      rrmx=rmx**2
      nrpt=0
      do 10 l=-mx,mx
      do 20 i=1,3
   20 v(i,1)=dble(l)*r(i,1)
      do 10 m=-mx,mx
      do 30 i=1,3
   30 v(i,2)=v(i,1)+dble(m)*r(i,2)
      do 10 n=-mx,mx
      do 40 i=1,3
   40 v(i,3)=v(i,2)+dble(n)*r(i,3)
      vv=1d10
      do 70 k=1,27
   70 vv=min(vv,
     &    (v(1,3)-p(1,k))**2+(v(2,3)-p(2,k))**2+(v(3,3)-p(3,k))**2)
      if(vv .lt. rrmx) then
      nrpt=nrpt+1
      if(nrpt .gt. nrptmx) call errtrp(1,'genlpt','nrpt too large')
      do 50 i=1,3
   50 rpt(i,nrpt)=v(i,3)
      endif
   10 continue
c     --- sort them in ascending order of the length
      call gpsort(rpt,nrpt)
      return
      end
