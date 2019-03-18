      function nfqlty(bzqlty,ibrav)
c-----------------------------------------------------------------------
c     Given the type of bravais lattice and the quality specification
c     This program returns appropriate nf parameter.
c     t: test run  l: low quality  m: medium  h: hifg  u: ultra high
c     also an integer number from 0 to 99 can be specified.
c
c     Now bzqlty can contains key-word 'fbz' in addition to the
c     above BZ quality specifications. If 'fbz' aooears, a full
c     Brillouin zone will be taken instead of an irreducible zone.
c
c     (0)unknown
c     (1)fcc (2)bcc (3)hcp (4)sc (5)bct (fct) (6)simple tetragonal
c     (7)face centered orthorhombic (8)body centered orthorhombic
c     (9)base centered orthorhonbic (10)simple orthorhombic
c     (11)base centered monoclinic (12)simple monoclinic
c     (13)triclinic (14)rhombohedral (trigonal) (15)fct (bct)
c     (16)fcc-like special definition 1 (sp1)
c     (17)bcc-like special definition 2 (sp2)
c             .....
c     (23)bco-like special definition 8 (sp8)
c
c     coded by H.Akai, April 1992,Osaka
c     revised by H.Akai, Feb. 1993, Osaka
c     modified by H. Akai, 11 Dec. 2014, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer nfdata(5,0:23)
      character qlist(5)*1,bzqlty*8,buff*8
      data nfdata/2,4,7,9,13,
     &            3,5,9,12,16, 3,6,11,15,20, 1,2,4,5,7, 2,4,7,9,13,
     &            2,4,8,10,15, 1,2,5,7,9, 0,0,0,0,0, 2,4,8,10,15,
     &            0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
     &            0,0,0,0,0, 3,5,10,12,17, 2,4,8,10,15,
     &            0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
     &            0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 2,4,7,9,11/
      data qlist/'t','l','m','h','u'/
      if(ibrav .lt. 0 .or. ibrav .gt. 23)
     &      call errtrp(1,'nfqlty','illegal ibrav')
      buff=bzqlty
      i=index(buff,'fbz')
      if(i .gt. 0) then
      buff(i:i+2)=' '
      m=lftfil(buff)
      endif
      do 10 i=1,5
      if(buff .eq. qlist(i)) then
      nfqlty=nfdata(i,ibrav)
      return
      endif
   10 continue
      read(buff,'(i4)')nfqlty
      if(nfqlty .ge. 0 .or. nfqlty .le. 1000) return
      call errtrp(1,'nfqlty','illegal bzqlty')
      end
