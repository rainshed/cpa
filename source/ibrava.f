      function ibrava(brvtyp)
c-----------------------------------------------------------------------
c     This program returns the index of the bravais lattice.
c     (1)fcc (2)bcc (3)hcp(hex) (4)sc (5)bct (fct) (6)simple tetragonal
c     (7)face centered orthorhombic (8)body centered orthorhombic
c     (9)base centered orthorhonbic (10)simple orthorhombic
c     (11)base centered monoclinic (12)simple monoclinic
c     (13)triclinic (14)rhombohedral (trigonal) (15)fct (bct)
c     (16)aux(prv) (primitive unit vector are to be read in)
c     for a practical reason, fct and bct are treated differently.
c     coded by H.Akai, April, 1992, Osaka
c     revised 26 Dec. 1994, Osaka
c     revised 26 Jan. 2017, Tokyo
c-----------------------------------------------------------------------
      character brvtyp*(*),brav(19)*3
      integer ibrav(19)
      data brav/'fcc','bcc','hcp','sc','bct','st','fco','bco',
     &          'bso','so','bsm','sm','trc','rhb','fct','trg','hex',
     &          'aux','prv'/,
     &     ibrav/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,14,3,16,16/
      ibrava=0
      do 10 i=1,19
   10 if(brvtyp .eq. brav(i)) ibrava=ibrav(i)
      end
