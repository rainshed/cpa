      subroutine madlng(eta,vc,atmicp,natm,rpt,nrpt,gpt,ngpt,rmt,amdlng
     &                 ,asa)
c-----------------------------------------------------------------------
c     Calculate constants needed for the constructin of the madelung
c     potential. Optimum value of 'eta' (in unit of 2*pi/a) for the
c     ewalt sum is used in this program.
c
c     eta:    convergence parameter
c     vc:     unit cell volume (in unit of a*3)
c     atmicp: position of atoms in the unit cell (in unit of a)
c     natm:   number of atoms in the unit cell
c     rpt:    neighbouring lattice points (in unit of a)
c     nrpt:   number of the lattice points given in 'rpt'
c     gpt:    reciprocal lattice vectors (in unit of 2*pi/a)
c     ngpt:   number of the reciprocal vectors given in 'gpt'
c     rmt:    muffin-tin radius of the i-th atom (in unit of a)
c     amdlng: calculated constants used to construct static potential
c
c     Elctro-static potential for i-th muffin-tin potential is
c     obtained by use of 'amdlng' as
c
c     v(i)=(contribution from inside the muffintin sphere)
c          +amdlng(i,1)*(q(1)/a)+amdlng(i,2)*(q(2)/a)+...
c          +amdlng(i,natm)*(q(natm)/a)
c
c     Here q(j) (j=1,natm) is the total charges, whichi include the
c     nuclear charge, inside the j-th muffin-tin sphere.
c     total static energy is also obtained as
c
c     w=(contribution from inside) - (1/2a) q(i)*amdlng(i,j)*q(j).
c
c     The sum over i and j should be taken. See H.A.'s note for
c     further details.
c
c     coded by H.Akai April 1992, Osaka
c     revised so as to fit also for ASA cases.
c     by H.Akai 14 Aug. 1995, Duisburg
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension atmicp(3,natm),rmt(natm),rpt(3,nrpt),gpt(3,ngpt)
     &         ,amdlng(natm,natm),dr(3)
      logical asa
      pi=4d0*atan(1d0)
      twopi=2d0*pi
c     --- optimal 'eta' parameter ---
      do 10 i=1,natm
      do 10 j=i,natm
      do 20 n=1,3
   20 dr(n)=atmicp(n,j)-atmicp(n,i)
c     --- reciprocal space sum ---
      sum=0d0
      do 30 k=1,ngpt
      arg=-twopi*(gpt(1,k)*dr(1)+gpt(2,k)*dr(2)+gpt(3,k)*dr(3))
      gg=gpt(1,k)**2+gpt(2,k)**2+gpt(3,k)**2
   30 sum=sum+cos(arg)*exp(-gg/4d0/eta**2)/gg
c     write(*,'(a,1p,e13.6)')
c    &  ' last reciprocal element=', cos(arg)*exp(-gg/4d0/eta**2)/gg
      amdlng(i,j)=-2d0*sum/vc/pi
c     --- real space sum ---
      sum=0d0
      do 40 k=1,nrpt
      arg=sqrt((rpt(1,k)+dr(1))**2+(rpt(2,k)+dr(2))**2
     &        +(rpt(3,k)+dr(3))**2)
   40 sum=sum+erfc(twopi*eta*arg)/arg
c     write(*,'(a,1p,e13.6)')
c    &  ' last real space element=', erfc(twopi*eta*arg)/arg
      amdlng(i,j)=amdlng(i,j)-2d0*sum
c     --- additional term coming from k=0 and that from r=0 ---
      if(i .eq. j) then
      amdlng(i,j)=amdlng(i,j)+8d0*sqrt(pi)*eta
      else
      arg=sqrt(dr(1)**2+dr(2)**2+dr(3)**2)
      amdlng(i,j)=amdlng(i,j)-2d0*erfc(twopi*eta*arg)/arg
      endif
      amdlng(i,j)=amdlng(i,j)+1d0/twopi/vc/eta**2
      amdlng(j,i)=amdlng(i,j)
   10 continue
      if(asa) return
c     write(6,1000)eta,amdlng(1,1)
c1000 format(' eta, b (see p.h.d. note) =',f12.5,f15.9)
c     --- interstitial volume ---
      sum=0d0
      do 50 i=1,natm
   50 sum=sum+rmt(i)**3
      volint=vc-4d0*pi*sum/3d0
      if(volint .lt. 0d0) call errtrp(1,'madlng','negative volint')
c     --- radius corresponding to the interstitial volume ---
      rintrs=(3d0*volint/4d0/pi)**(1d0/3d0)
      vov=volint/vc
c     --- now calculate the constant used for the static potential ---
      do 70 j=1,natm
      vjvorj=(rmt(j)/rintrs)**3/rmt(j)
      sum=0d0
      do 60 i=1,natm
      vivo=(rmt(i)/rintrs)**3
      vivor=vivo*vov/rmt(i)
      sum=sum+vivo*(-(3d0/5d0)*vivor+amdlng(i,j))
   60 amdlng(i,j)=amdlng(i,j)-3d0*(vivor+vjvorj)
      do 70 i=1,natm
   70 amdlng(i,j)=amdlng(i,j)+sum
c     ---finally i convert amdlng into another symmetric matrix which
c        is more useful. see h.a. note.
      do 80 i=1,natm
      sum=0d0
      do 90 k=1,natm
   90 sum=sum+amdlng(i,k)*(rmt(k)/rintrs)**3
      do 80 j=1,natm
   80 amdlng(i,j)=amdlng(i,j)+sum
      return
      end
