      subroutine genrpt(r,rpt,nrptmx,nrpt,rmx,p,natm)
c----------------------------------------------------------------------
c     Given three primitive translation vectors, this program generates
c     a set of lattice vectors, whose maximum length do not exceeds
c     'rmx'. The vectors are stored in rpt in ascending order of the
c     length.
c     coded by H.Akai April 1992, Osaka
c     tiny modification by H.Akai, 6 Dec. 1995, Osaka
c     major modification by H.Akai, 1 Jul. 2007, Osaka
c     Modified so that the set of R-points preserves the crystal
c     symmetry.
c     structure.
c     H. Akai, 30 Nov. 2014, Tokyo
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mx=50)
      dimension r(3,3),rpt(3,nrptmx),v(3,3),p(3,natm)
c   
c     Generate the lattice points for which the distances from all
c     the atomic positions in the primitive cell are less than rmx.
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
      do 70 k1=1,natm
      do 70 k2=1,natm
      rr=0d0
      do 80 i=1,3
   80 rr=rr+(v(i,3)-p(i,k1)+p(i,k2))**2
   70 vv=min(vv,rr)
      if(vv .lt. rrmx) then
      nrpt=nrpt+1
      if(nrpt .gt. nrptmx) call errtrp(1,'genrpt','nrpt too large')
      do 50 i=1,3
   50 rpt(i,nrpt)=v(i,3)
      endif
   10 continue
c     --- sort them in ascending order of the length
      call gpsort(rpt,nrpt)
      end
