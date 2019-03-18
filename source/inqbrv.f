      function inqbrv(anclr)
c----------------------------------------------------------------------
c     Given the nuclear charges this function returns guess type of the
c     bravais lattice after experimental data.
c     'istrct' is taken from sargent-welch table. 1..8 correspond
c     to (1)fcc (2)bcc (3)sc (4)hcp (5)rhb (6)st (7)so (8)sm
c     'ibra' establishes the link between sargent-welch notification
c     and 'ibrav' used in the KKR program.
c     coded by H.Akai, 25 Aug. 1991.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension istrct(100),ibra(8)
      data istrct/4,4,2,4,5,4,4,3,3,1,2,4,1,1,8,7,7,1,2,1,4,4,
     &            2,2,1,2,4,1,1,4,7,1,5,4,7,1,2,1,4,4,2,2,4,4,
     &            1,1,1,4,6,6,5,4,7,1,2,2,4,1,4,4,4,5,2,4,4,4,
     &            4,4,4,1,4,4,2,2,4,4,1,1,1,5,4,1,5,8,0,1,2,2,
     &            1,1,7,7,7,8,4,0,0,0,0,0/
      data ibra/1,2,4,3,14,6,10,12/
      nclr=anclr
      if(nclr .lt. 1 .or. nclr .gt. 100) then
      write(6,1000)
 1000 format('   ***err in guessb...anclr is illegal')
      stop
      endif
      inqbrv=ibra(istrct(nclr))
      return
      end
