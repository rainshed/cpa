      subroutine brvsym(ibrav,isymop,ls,r,g,protat)
c-----------------------------------------------------------------------
c     Given type of the Bravais lattice (1,..,14 and 15), this program
c     returns the compatible point symmetry operations. Correspondence
c     between ibrav and the type of the Bravais lattice is given in the
c     subroutines 'prmvec', 'ibrava' and 'bravai'.
c     Coded by H.Akai, 14 Jan. 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 r(3,3),g(3,3),protat(9,48),s(3,3)
      integer isym(24,7,2),isymop(24,2),isystm(17)
c--- without spin-orbit coupling (ls=0)
      data isym/
c        1:cubic  2:hexagonal  3:tetragonal
     &   24*1,    12*1,12*0,   1,0,0,4*1,0,0,1,8*0,1,1,4*0,
c        4:monoclinic   5:triclinic   6:orthorhombic
     &   1,4*0,1,18*0,  1,23*0,       1,3*0,3*1,17*0,
c        7:trigonal
     &   1,0,1,0,1,0,1,0,1,0,1,0,12*0,
c--- with spin-orbit coupling (ls=1)
c        1:cubic
     &   1,0,0,1,0,0,1,0,0,1,8*0,0,0,4*0,
c        2:hexagonal
     &   6*1,6*0,12*0,     
c        3:tetragonal
     &   1,0,0,1,0,0,1,0,0,1,8*0,0,0,4*0,
c        4:monoclinic   5:triclinic   6:orthorhombic
     &   1,4*0,1,18*0,  1,23*0,       1,3*0,0,0,1,17*0,
c        7:trigonal
     &   1,0,1,0,1,0,0,0,0,0,0,0,12*0/
c
      data isystm
c        system of each Bravais lattice
     &  /1,1,2,1,3,3,6,6,6,6,4,4,5,7,3,1,2/
      data zero/1d-6/
      if(ls .eq. 0) then
c     --- no spin-orbit coupling assumed
      ind=1
      else if(ls .eq. 1) then
c     --- spin-orbit coupling assumed
      ind=2
      else
      call errtrp(1,'brvsym','illegal ls')
      endif
      do 10 i=1,24
   10 isymop(i,1)=isym(i,isystm(ibrav),ind)
      if(ibrav .eq. 16 .or. ibrav .eq. 17) then
c     --- Bravai lattice is not exactly known
c     --- ibrav=16 means neither hexagonal nor trigonal
c     --- ibrav=17 means either hexagonal or trigonal
c     write(*,'(1x,9f8.5)')r
      do 20 iop=1,24
      if(isymop(iop,1) .ne. 0) then
c     --- check if this rotataion maps the Bravai lattice to itself.
c     --- rotate the cartesian cordinates of primitive vectors.
c     --- r -> s
      call rotatm(r,s,protat(1,iop),3)
c     write(*,'(1x,9f8.5)')s
      do 30 i=1,3
      do 30 j=1,3
c     --- s is expressed as a sum of the primitive vectors r1, r2
c         and r3 as s=a1*r1+a2*r2+a3*r3. the inner product s*g1, s*g2
c         and s*g3 give the coeficients a1, a2 and a3.
      a=s(1,i)*g(1,j)+s(2,i)*g(2,j)+s(3,i)*g(3,j)
c     write(*,'(1x,e15.7)')a
c     --- check if a is an integer.
      b=dble(int(a+1000.5d0))-1d3
c     --- b is the nearest integer of a.
      if(abs(a-b) .gt. zero) then
c     write(*,*)a,b
c     --- s cannot be expressed as any combinations of r1, r2 and r3.
      isymop(iop,1)=0
      go to 20
      endif
   30 continue
c     --- all s are expressed as a combination of primitive vectors r,
c         meaning that r are mapped to itself by this rotation.
      endif
   20 continue
      endif
      do 40 i=1,24
   40 isymop(i,2)=isymop(i,1)
      end
