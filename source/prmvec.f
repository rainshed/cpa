      subroutine prmvec(ibrav,coa,boa,alpha,beta,gamma,vc,r,g,angl)
c-----------------------------------------------------------------------
c     Given the type of the bravais lattice 'ibrav', c/a, b/a ratio,
c     and alpha, beta and gamma, this program generates the primitive
c     translation vectors r(i,j) and those of the reciprocal lattice
c     g(i,j) (i specifies x, y and z-components, j the independent
c     vectors), and the primitive cell volume, 'vc'. Specification
c     by 'ibrav' is preferentially activated irrespective of the
c     values of coa, boa, etc. 'fct' is treated as an independent
c     structure. Thus, 15 type of lattice appear instead of 14.
c
c     (1)fcc (2)bcc (3)hcp (4)sc (5)bct (fct) (6)simple tetragonal
c     (7)face centered orthorhombic (8)body centered orthorhombic
c     (9)base centered orthorhonbic (10)simple orthorhombic
c     (11)base centered monoclinic (12)simple monoclinic
c     (13)triclinic (14)rhombohedral (trigonal) (15)fct (bct)
c     (16)aux (primitive vectors are to be read in)
c     (17)prv (primitive vectors are to be read in)

c
c     coded by H.Akai, April, 1992, Osaka
c     revised 26 Dec. 1994, Osaka
c     revised 12 Jan. 1997, Osaka
c     triclinic revised by K. Yamamoto, 1999/02/1
c     monoclinic and triclinic revised by M. Ogura, Apr. 2013, Osaka
c
c     Modified so as to follow obverse rhombohedral coordinates
c     instead of the reverse setting used in the previous versions.
c     by H. Akai, 24 Dec. 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(c0=0d0,c1=0.5d0,c2=-0.5d0,c3=1d0,c4=-1d0,
     &          c5=8.660254037844386d-1,c6=-c5)
c     --- c5=sqrt(3d0)/2d0 ---
      real*8 r(3,3),g(3,3),prmr(3,3,15),rlng(3),rprd(3),angl(3),rtilt(3)
      data prmr
     &  /c0,c1,c1, c1,c0,c1, c1,c1,c0,  c2,c1,c1, c1,c2,c1, c1,c1,c2,
     &   c1,c6,c0, c1,c5,c0, c0,c0,c3,  c3,c0,c0, c0,c3,c0, c0,c0,c3,
     &   c2,c1,c1, c1,c2,c1, c1,c1,c2,  c3,c0,c0, c0,c3,c0, c0,c0,c3,
     &   c0,c1,c1, c1,c0,c1, c1,c1,c0,  c2,c1,c1, c1,c2,c1, c1,c1,c2,
     &   c1,c2,c0, c1,c1,c0, c0,c0,c3,  c3,c0,c0, c0,c3,c0, c0,c0,c3,
     &   c1,c2,c0, c1,c1,c0, c0,c0,c3,  c3,c0,c0, c0,c3,c0, c0,c0,c3,
     &   c3,c0,c0, c0,c3,c0, c0,c0,c3,  c0,c0,c0, c0,c0,c0, c0,c0,c0,
     &   c0,c1,c1, c1,c0,c1, c1,c1,c0/
      data zero/1d-6/
c
      pi=4d0*atan(1d0)
      rad=pi/180d0
c     rad=0d0
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,160),
     &       ibrav
      call errtrp(1,'prmvec','illegal ibrav')
c     --- fcc case ---
   10 boa=1d0
      coa=1d0
      alpha=90d0
      beta=90d0
      gamma=90d0
      go to 170
c     --- bcc case ---
   20 boa=1d0
      coa=1d0
      alpha=90d0
      beta=90d0
      gamma=90d0
      go to 170
c     --- hcp case ---
   30 boa=1d0
      if(abs(coa) .lt. 1d-6) coa=sqrt(8d0/3d0)
      alpha=90d0
      beta=90d0
      gamma=120d0
      go to 170
c     --- sc case ---
   40 boa=1d0
      coa=1d0
      alpha=90d0
      beta=90d0
      gamma=90d0
      go to 170
c     --- bct (fct) case ---
   50 boa=1d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      alpha=90d0
      beta=90d0
      gamma=90d0
      go to 170
c     --- st case ---
   60 boa=1d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      alpha=90d0
      beta=90d0
      gamma=90d0
      go to 170
c     --- fco case ---
   70 alpha=90d0
      beta=90d0
      gamma=90d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      if(abs(boa) .lt. 1d-6) boa=1d0
      go to 170
c     --- bco case ---
   80 alpha=90d0
      beta=90d0
      gamma=90d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      if(abs(boa) .lt. 1d-6) boa=1d0
      go to 170
c     --- bso ---
   90 alpha=90d0
      beta=90d0
      gamma=90d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      if(abs(boa) .lt. 1d-6) boa=1d0
      go to 170
c     --- so case ---
  100 alpha=90d0
      beta=90d0
      gamma=90d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      if(abs(boa) .lt. 1d-6) boa=1d0
      go to 170
c     --- bsm case ---
  110 continue
      alpha=90d0
      gamma=90d0
      if(abs(beta) .lt. 1d-6) beta=90d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      if(abs(boa) .lt. 1d-6) boa=1d0
      r(1,1)=5d-1
      r(2,1)=-5d-1*boa
      r(3,1)=0d0
      r(1,2)=5d-1
      r(2,2)=5d-1*boa
      r(3,2)=0d0
      r(1,3)=cos(pi/180d0*beta)*coa
      r(2,3)=0d0
      r(3,3)=sqrt(1d0-(cos(pi/180d0*beta))**2)*coa
      go to 180
c     --- sm case ---
  120 continue
      alpha=90d0
      gamma=90d0
      if(abs(beta) .lt. 1d-6) beta=90d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      if(abs(boa) .lt. 1d-6) boa=1d0
      r(1,1)=1d0
      r(2,1)=0d0
      r(3,1)=0d0
      r(1,2)=0d0
      r(2,2)=boa
      r(3,2)=0d0
      r(1,3)=cos(pi/180d0*beta)*coa
      r(2,3)=0d0
      r(3,3)=sqrt(1d0-(cos(pi/180d0*beta))**2)*coa
      go to 180
c     --- trc case ---
  130 continue
      if(abs(alpha) .lt. 1d-6) alpha=90d0
      if(abs(beta) .lt. 1d-6) beta=90d0
      if(abs(gamma) .lt. 1d-6) gamma=90d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      if(abs(boa) .lt. 1d-6) boa=1d0
      cw=1d0-(cos(pi/180d0*alpha))**2-(cos(pi/180d0*beta))**2-
     &   (cos(pi/180d0*gamma))**2+2*cos(pi/180d0*alpha)*
     &   cos(pi/180d0*beta)*cos(pi/180d0*gamma)
      r(1,1)=1d0
      r(2,1)=0d0
      r(3,1)=0d0
      r(1,2)=cos(pi/180d0*gamma)*boa
      r(2,2)=sin(pi/180d0*gamma)*boa
      r(3,2)=0d0
      r(1,3)=cos(pi/180d0*beta)*coa
      r(2,3)=(cos(pi/180d0*alpha)-cos(pi/180d0*beta)*
     &   cos(pi/180d0*gamma))/sqrt(1d0-(cos(pi/180d0*gamma))**2)*coa
      r(3,3)=sqrt(cw)/sqrt(1d0-(cos(pi/180d0*gamma))**2)*coa
      go to 180
c     --- rhb (trg) case ---
c     --- obverse setting is used
c         Note: in the previous version, the reverse setting is used.
  140 boa=1d0
      coa=1d0
      if(abs(alpha) .lt. 1d-6) alpha=60d0
      beta=alpha
      gamma=alpha
      theta=2d0*pi*alpha/360d0
      z=sqrt((5d-1+cos(theta))/(1d0-cos(theta)))
      snrm=sqrt(1d0+z**2)
      r(1,1)=(sqrt(3d0)/2d0)/snrm
      r(2,1)=-(1d0/2d0)/snrm
      r(3,1)=z/snrm
      r(1,2)=0d0
      r(2,2)=1d0/snrm
      r(3,2)=z/snrm
      r(1,3)=-(sqrt(3d0)/2d0)/snrm
      r(2,3)=-(1d0/2d0)/snrm
      r(3,3)=z/snrm
c     ---following gives a different choice of the axis orientation.
c        In this case one of the axis directs to [100].
c     theta=4d0*atan(1d0)*alpha/180d0
c     r(1,1)=1d0
c     r(2,1)=0d0
c     r(3,1)=0d0
c     r(1,2)=cos(theta)
c     r(2,2)=sin(theta)
c     r(3,2)=0d0
c     r(1,3)=cos(theta)
c     r(2,3)=cos(theta)*tan(theta/2d0)
c     r(3,3)=sqrt(1d0-r(1,3)**2-r(2,3)**2)
      go to 180
c     --- fct (bct) case ---
  150 boa=1d0
      if(abs(coa) .lt. 1d-6) coa=1d0
      alpha=90d0
      beta=90d0
      gamma=90d0
      go to 170
  160 do 190 i=1,3
  190 rlng(i)=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
      do 200 i=1,3
      j=mod(i,3)+1
  200 rprd(i)=(r(1,i)*r(1,j)+r(2,i)*r(2,j)+r(3,i)*r(3,j))
     &        /rlng(i)/rlng(j)
      if(    abs(rlng(1)-rlng(2)) .lt. zero
     & .and. abs(rlng(2)-rlng(3)) .lt. zero
c     --- a=b=c
     & .and. abs(rprd(1)-rprd(2)) .lt. zero
     & .and. abs(rprd(2)-rprd(3)) .lt. zero
c     --- alpha=beta=gamma
     & .and. abs(r(3,1)-r(3,2)) .lt. zero
     & .and. abs(r(3,2)-r(3,3)) .lt. zero) then
c     --- z-component of a, b, c is zero
c     --- this case should be trigonal ---
      ibrav=17
      coa=1d0
      boa=1d0
c     write(*,'(1x,a)')'this case is trigonal'
      go to 180
      endif
      if(    abs(rlng(1)-rlng(2)) .lt. zero
c     --- a=b
     &     .and. (abs(rprd(1)-5d-1) .lt. zero
c     --- gamma=60
     &     .or. abs(rprd(1)+5d-1) .lt .zero)
c     --- gamma=120
     &     .and. abs(rprd(2)) .lt. zero
     &     .and. abs(rprd(3)) .lt. zero) then
c     --- alpha=beta=90
c     --- this case should be hcp ---
c     write(*,'(1x,a)')'this case is hexagonal'
      ibrav=17
      coa=2d0*r(3,3)
      boa=1d0
      go to 180
      endif
c     --- at the moment not fixed but likeky neither hex nor trg---
      ibrav=16
      go to 180
  170 continue
      do 210 i=1,3
      r(1,i)=prmr(1,i,ibrav)
      r(2,i)=prmr(2,i,ibrav)*boa
  210 r(3,i)=prmr(3,i,ibrav)*coa
c     --- calculate primitive cell volume ---
  180 do 220 i=1,3
      call vrotat(rtilt,r(1,i),rad*angl(3),rad*angl(2),rad*angl(1))
      do 220 j=1,3
  220 r(j,i)=rtilt(j)
      vc=(r(2,1)*r(3,2)-r(3,1)*r(2,2))*r(1,3)
     &  +(r(3,1)*r(1,2)-r(1,1)*r(3,2))*r(2,3)
     &  +(r(1,1)*r(2,2)-r(2,1)*r(1,2))*r(3,3)
c     --- primitive reciprocal lattice vectors ---
      g(1,1)=(r(2,2)*r(3,3)-r(3,2)*r(2,3))/vc
      g(2,1)=(r(3,2)*r(1,3)-r(1,2)*r(3,3))/vc
      g(3,1)=(r(1,2)*r(2,3)-r(2,2)*r(1,3))/vc
      g(1,2)=(r(2,3)*r(3,1)-r(3,3)*r(2,1))/vc
      g(2,2)=(r(3,3)*r(1,1)-r(1,3)*r(3,1))/vc
      g(3,2)=(r(1,3)*r(2,1)-r(2,3)*r(1,1))/vc
      g(1,3)=(r(2,1)*r(3,2)-r(3,1)*r(2,2))/vc
      g(2,3)=(r(3,1)*r(1,2)-r(1,1)*r(3,2))/vc
      g(3,3)=(r(1,1)*r(2,2)-r(2,1)*r(1,2))/vc
      vc=abs(vc)
      end
