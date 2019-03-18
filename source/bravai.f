      function bravai(i)
c-----------------------------------------------------------------------
c     This program returns the name of the bravais lattice.
c     (1)fcc (2)bcc (3)hcp (4)sc (5)bct (fct) (6)simple tetragonal
c     (7)face centered orthorhombic (8)body centered orthorhombic
c     (9)base centered orthorhonbic (10)simple orthorhombic
c     (11)base centered monoclinic (12)simple monoclinic
c     (13)triclinic (14)rhombohedral (trigonal) (15)fct (bct)
c     (16)primitive vectors are to be read in (aux)
c     (17)primitive vectors are to be read in (aux)
c     for a practical reason, fct and bct are treated differently.
c     coded by H.Akai, April, 1992, Osaka
c-----------------------------------------------------------------------
      character*3 bravai,brav(17)
      data brav/'fcc','bcc','hcp','sc','bct','st','fco','bco',
     &          'bso','so','bsm','sm','trc','rhb','fct',
     &          'aux','aux'/
      if(i .lt. 1 .or. i .gt. 17) then
      bravai=' '
      return
      else
      bravai=brav(i)
      return
      endif
      end
