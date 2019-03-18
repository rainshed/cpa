      subroutine corcnf(anclr,config,corlvl,ncmpx,ef)
c----------------------------------------------------------------------
c     Given a nuclear charge, this program returns the core state
c     configuration 'config', 18*2 array indicating the occupation
c     number of the states specifies by (n,l) and the spin.
c
c      j    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
c     -----------------------------------------------------------
c      nl  1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
c
c     coded by H.Akai, 28 July 1993, Osaka.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 anclr(ncmpx),config(18,ncmpx,2),corlvl(18,ncmpx,2),oc(4)
     &      ,ef(2)
      integer ndd(0:107),nconf(15),ncrit(6),nc(6),nval(15,2)
c
c     ndd gives the number of the core shells that are not included in 
c     a rare-gas core state. - sign denotes that the core state should
c     be treated as a semi-cores, i.e. they should not be treated as a
c     core state but a part of the valence states..
c
      data ndd
c       ---------------------------------------------------------------
c        Vc     H      He     Li     Be     B      C      N      O
     & / 0,     0,     0,     0,     0,     0,     1,     1,     1,
c       ---------------------------------------------------------------
c        F      Ne     Na     Mg     Al     Si     P      S      Cl
     &   1,     0,     0,     0,     0,     0,     0,     0,     1,
c       ---------------------------------------------------------------
c        Ar     K      Ca     Sc     Ti     V      Cr     Mn     Fe
     &   0,     0,     0,     0,     0,     0,     0,     0,     0,
c       ---------------------------------------------------------------
c        Co     Ni     Cu     Zn     Ga     Ge     As     Se     Br
     &   0,     0,     0,     0,     1,     1,     2,     2,     2,
c       ---------------------------------------------------------------
c        Kr     Rb     Sr     Y      Zr     Nb     Mo     Tc     Ru
     &   0,     0,     0,     0,     0,     0,     0,     0,     0,
c       ---------------------------------------------------------------
c        Rh     Pd     Ag     Cd     In     Sn     Sb     Te     I
     &   0,     0,     0,     0,     1,     1,     2,     2,     1,
c       ---------------------------------------------------------------
c        Xe     Cs     Ba     La     Ce     Pr     Nd     Pm     Sm
     &   0,     0,     0,     0,     1,     1,     1,     1,     1,
c       ---------------------------------------------------------------
c        Eu     Gd     Tb     Dy     Ho     Er     Tm     Yb     Lu
     &   1,     1,     1,     1,     1,     1,     1,     1,     1,
c       ---------------------------------------------------------------
c        Hf     Ta     W      Re     Os     Ir     Pt     Au     Hg
     &   1,     1,     1,     1,     1,     1,     1,     1,     1,
c       ---------------------------------------------------------------
c        Tl     Pb     Bi     Po     At     Rn     Fr     Ra     Ac
     &   2,     2,     3,     3,     3,     0,     0,     0,     0,
c       ---------------------------------------------------------------
c        Th     Pa     U      Np     Pu     Am     Cm     Bk     Cf
     &   0,     0,     0,     0,     0,     0,     0,     0,     0,
c       ---------------------------------------------------------------
c        Es     Fm     Md     No     Lr     Unq    Unp    Unh    Uns
     &   0,     0,     0,     0,     0,     1,     1,     1,     1   /
c       ---------------------------------------------------------------
c
      data nconf/1,1,3,1,3,5,1,3,5,1,3,7,5,1,3/, oc/7d0,5d0,1d0,3d0/
     &    ,ncrit/2,10,18,36,54,86/, nc/1,3,5,8,11,15/
c
      data nval
c         -----------------------------------------------
c          La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
     &    / 3, 4, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3,
c    &    / 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3,
c         -----------------------------------------------
c          Ac Th Pa  U Np Pu Am Cm Bk Cf Es Fm Md No Lr
     &      3, 4, 3, 3, 3, 2, 2, 3, 2, 2, 2, 2, 2, 2, 3/
c         -----------------------------------------------
c
      do 60 j=1,ncmpx
      nz=nint(anclr(j))
      nd=0
      do 10 i=1,6
   10 if(nz .ge. ncrit(i)) nd=nc(i)+ndd(nz)
      do 20 i=1,nd
   20 config(i,j,1)=dble(nconf(i))
      do 40 i=nd+1,18
   40 config(i,j,1)=0d0
      do 50 i=1,18
   50 config(i,j,2)=config(i,j,1)
c
c     ---lanthanide
      if(nz .ge. 57 .and. nz. le. 71 .and. config(12,j,1) .gt. 0d0) then
      nf=nz-54-nval(nz-56,1)
c     --- Be careful about the direction chosen for the f spin
      dv=0d0
      do 70 i=1,11
   70 dv=dv+dble(nconf(i))*(ef(1)-corlvl(i,j,1)-ef(2)+corlvl(i,j,2))
      if(dv .gt. 0d0) then
c     write(*,*)'case 1',ef(1)+corlvl(1,j,1),ef(2)+corlvl(1,j,2)
      config(12,j,1)=dble(min(7,nf))
      config(12,j,2)=dble(max(0,nf-7))
c     config(12,j,1)=7d0
c     config(12,j,2)=6.8d0
      else
c     --- up/down is reversed ---
c     write(*,*)'case 2',ef(1)+corlvl(1,j,1),ef(2)+corlvl(1,j,2)
      config(12,j,2)=dble(min(7,nf))
      config(12,j,1)=dble(max(0,nf-7))
c     config(12,j,2)=7d0
c     config(12,j,1)=6.8d0
      endif
      endif
c
c     ---actinide
      if(nz .ge. 89 .and. nz. le. 103 .and. 
     &              config(16,j,1) .gt. 0d0) then
      nf=nz-86-nval(nz-88,2)
c     --- Be careful about the direction chosen for the f spin
      dv=0d0
      do 80 i=1,15
   80 dv=dv+dble(nconf(i))*(ef(1)-corlvl(i,j,1)-ef(2)+corlvl(i,j,2))
      if(dv .gt. 0d0) then
      config(16,j,1)=dble(min(7,nf))
      config(16,j,2)=dble(max(0,nf-7))
      else
c     --- up/down is reversed ---
      config(16,j,2)=dble(min(7,nf))
      config(16,j,1)=dble(max(0,nf-7))
      endif
      endif
      write(*,'(a,i3)') '   core configuration for Z=',nz
      write(*,'(2a)') '   state  ',
     & '1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s'
      write(*,'(a,18i3)')'    up   ',(int(config(i,j,1)),i=1,18)
   60 write(*,'(a,18i3/)')'   down  ',(int(config(i,j,2)),i=1,18)
c     stop
      end
