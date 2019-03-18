      subroutine mkemap(id,atmicp,itype,r,natm,dlim,dist
     &                 ,ityp,ideg,nb,nd)
c-----------------------------------------------------------------------
c     Make a map of near neighbors. It gives types and the distances of
c     the neighboring atoms for each type of atoms.
c     coded by H.Akai April, 1992
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 atmicp(3,natm),r(3,3),dist(nb),sft(3)
      integer itype(natm),ityp(nb),ideg(nb)
      data small/1d-6/ , nc/3/
c     --- fix fc such that workarea will not overflow
c     nc=5d-1*((dble(nb)/dble(natm))**(1d0/3d0)-1d0)+1d-6
c     --- first fill in null data.
      do 100 j=1,nb
      dist(j)=0d0
      ityp(j)=0
  100 ideg(j)=0
c     --- look for an atom with specified type.
      do 10 j=1,natm
      ifit=j
      if(itype(j) .eq. id) go to 20
   10 continue
      call errtrp(2,'mkemap','type not appears')
      write(6,'(a,i2)')'   type=',id
      ityp(1)=id
      nd=0
      return
c
   20 ind=0
      do 30 i1=-nc,nc
      do 30 i2=-nc,nc
      do 30 i3=-nc,nc
      do 40 l=1,3
   40 sft(l)=dble(i1)*r(l,1)+dble(i2)*r(l,2)+dble(i3)*r(l,3)
      do 30 j=1,natm
      d=0d0
      do 50 l=1,3
   50 d=d+(atmicp(l,j)+sft(l)-atmicp(l,ifit))**2
      if(d .lt. dlim**2) then
      ind=ind+1
      if(ind .gt. nb) call errtrp(1,'mkemap','table overflows')
      dist(ind)=sqrt(d)
      ityp(ind)=itype(j)
      endif
   30 continue
c     write(6,1200)ind
c1200 format(' number of neighbors considered=',i4)
c     --- sort in accending order of both dist and ityp.
      do 60 j=2,ind
      a=dist(j)
      ia=ityp(j)
      do 70 i=j-1,1,-1
      if(dist(i) .lt. a) go to 80
      if(abs(dist(i)-a) .lt. small .and. ityp(i) .lt. ia) go to 80
      dist(i+1)=dist(i)
      ityp(i+1)=ityp(i)
   70 continue
      i=0
   80 dist(i+1)=a
      ityp(i+1)=ia
   60 continue
c     --- group the euivalent neighbors.
      nd=1
      ideg(1)=1
      do 90 j=2,ind
      if(abs(dist(j)-dist(j-1)) .gt. small
     &           .or. ityp(j) .ne. ityp(j-1)) then
      nd=nd+1
      dist(nd)=dist(j)
      ityp(nd)=ityp(j)
      endif
   90 ideg(nd)=ideg(nd)+1
c     write(6,2000)(dist(i),ityp(i),ideg(i),i=1,nd)
c2000 format((' dist=',f12.5,' typ=',i1,' deg=',i2))
      return
      end
