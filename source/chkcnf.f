      subroutine chkcnf(anclr,config,ncmpx,tof,mmxl,mxlcmp)
c----------------------------------------------------------------------
c     Given a nuclear charge, core-configuration and core energy level
c     this program checks if the core-configuration is all right in the
c     case of an open f-shell structure is used. The criterion of
c     a sensible core-configuration is that the spin moment of the
c     core f-states is parallel to the local spin moment of the atom.
c     If this is not the case, spin up and down are reversed in the
c     core-configurations.
c
c      j    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
c     -----------------------------------------------------------
c      nl  1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
c
c     coded by H.Akai, 28 Sep. 2013, Tokyo.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 anclr(ncmpx),config(18,ncmpx,2),tof(mmxl,ncmpx,2)
      integer mxlcmp(ncmpx)
c
      do 10 j=1,ncmpx
      nz=nint(anclr(j))
c     ---lanthanide
      if(nz .ge. 57 .and. nz. le. 71) then
      s=0d0
      do 20 is=1,2
      do 20 l=1,min(9,mxlcmp(j)**2)
   20 s=s+tof(l,j,is)*(-1d0)**(is-1)
      if(s*(abs(config(12,j,1))-abs(config(12,j,2))) .lt. 0d0)
     &    call swparr(config(12,j,1),config(12,j,2),1)
      endif
c     ---actinide
      if(nz .ge. 89 .and. nz. le. 103) then
      do 30 is=1,2
      do 30 l=1,min(9,mxlcmp(j)**2)
   30 s=s+tof(l,j,is)*(-1d0)**(is-1)
      if(s*(abs(config(16,j,1))-abs(config(16,j,2))) .lt. 0d0)
     &    call swparr(config(16,j,1),config(16,j,2),1)
      endif
c     write(*,'(a,i3)') '   core configuration for Z=',nz
c     write(*,'(2a)') '   state  ',
c    & '1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s'
c     write(*,'(a,18i3)')'    up   ',(int(config(i,j,1)),i=1,18)
c  10 write(*,'(a,18i3/)')'   down  ',(int(config(i,j,2)),i=1,18)
   10 continue
c     stop
      end
