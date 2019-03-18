      subroutine brvsy2(isys,isymop,ls,r,g,protat)
c-----------------------------------------------------------------------
c     New version of brvsym. The old version does not work always
c     correctly when the crystal is tilted. In this version the symmetry
c     of bravais lattice itself is checked without relying on the
c     information given by ibrav. All the symmetry operation is tried
c     when spin-orbit coupling is suppressed. When spin-orbit coupling
c     exists, only the rotations around z-axis are allowed.
c     isym contains the following data:
c
c     Without spin-orbit coupling:
c        cubic: 24*1,  hexagonal: 12*1, 12*0
c     With spin-orbit couping:
c        cubic: 1,0,0,1,0,0,1,0,0,1,14*0,  hegonal:6*1,18*0
c
c     Modified by H. Akai, Tokyo, Dec. 11, 2017 from
c     old version of brvsym coded by H.Akai, 14 Jan. 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 r(3,3),g(3,3),protat(9,48),s(3,3)
      integer isym(24,2,2),isymop(24,2)
c--- without spin-orbit coupling (ls=0)
      data isym/24*1,12*1,12*0,1,0,0,1,0,0,1,0,0,1,14*0,6*1,18*0/
      data zero/1d-6/
      if(ls .eq. 0) then
c     --- no spin-orbit coupling assumed
      ind=1
      else if(ls .eq. 1) then
c     --- spin-orbit coupling assumed
      ind=2
      else
      call errtrp(1,'brvsy2','illegal ls')
      endif
      do 10 i=1,24
   10 isymop(i,1)=isym(i,isys,ind)
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
      do 40 i=1,24
   40 isymop(i,2)=isymop(i,1)
      end
