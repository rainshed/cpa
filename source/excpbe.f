      subroutine excpbe(ro,v,w,dr,xr,msr,mna,meshr)
c-----------------------------------------------------------------------
c     driver routine for pbe gga subroutines
c     Adapted to kkr-code by H. Ebert
c     Modified by H. Akai, Nov. 2007
c     Last modified by H. Akai, Mar. 2013
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      real*8 ro(msr,mna,2),v(msr,mna,2),w(msr,mna,*)
     &      ,dr(msr,mna),xr(msr,mna)
      pi=4d0*atan(1d0)
      thrd=1d0/3d0
      thrd2=2d0/3d0
      conf = (3d0*pi**2)**thrd
      conrs = (3d0/(4d0*pi))**thrd
      fourpi=4d0*pi
      do 10 ia=1,mna
      do 10 k = 1,meshr
      exc = 0d0
      vxcup = 0d0
      vxcdn = 0d0
      if(ro(k,ia,1) .gt. 1d-20 .and.
     &       .not. abs(ro(k,ia,2)) .gt. ro(k,ia,1)) then
c     --- preparation
      rho=ro(k,ia,1)/fourpi
      sho=ro(k,ia,2)/fourpi
      if(k .le. meshr-1) then
      call drivtv(ro(1,ia,1),meshr-1,xr(1,ia),dr(1,ia),k,rho1,rho2)
      call drivtv(ro(1,ia,2),meshr-1,xr(1,ia),dr(1,ia),k,sho1,sho2)
      else
      rho1=0d0
      rho2=0d0
      sho1=0d0
      sho2=0d0
      endif
      rho1=rho1/fourpi
      rho2=rho2/fourpi
      sho1=sho1/fourpi
      sho2=sho2/fourpi
      rhou=5d-1*(rho+sho)
      rhod=5d-1*(rho-sho)
      rhou1=5d-1*(rho1+sho1)
      rhod1=5d-1*(rho1-sho1)
      rhou2=5d-1*(rho2+sho2)
      rhod2=5d-1*(rho2-sho2)
c --- begin the spin loop for exchange
      do 20 isp = 1,2
      if (isp .eq. 1) then
      d = 2d0*rhou
      dp = 2d0*rhou1
      dpp = 2d0*rhou2
      else
      d = 2d0*rhod
      dp = 2d0*rhod1
      dpp = 2d0*rhod2
      end if
      if(d .gt. 1d-20) then
      fk = conf*d**thrd
      ss = abs(dp)/(d*2d0*fk)
      uu = dp*(-dpp)/(d**2*(2d0*fk)**3)
      vv = (2d0*dp/xr(k,ia)+dpp)/(d*(2d0*fk)**2)
      call excpbex(d,ss,uu,vv,ex,vx)
      exc = exc + ex*(d/2d0)/rho
      else
      vx=0d0
      endif
      if (isp .eq. 1) vxcup = vx
   20 if (isp .eq. 2) vxcdn = vx
c --------------------------------------------------------- correlation
      d = rho
      zet = (rhou-rhod)/rho
      rs = conrs/d**thrd
      fk = 1.91915829d0/rs
      sk = dsqrt(4d0*fk/pi)
      g = ((1d0+zet)**thrd2+(1d0-zet)**thrd2)/2d0
      tt = abs(rho1)/(d*2d0*sk*g)
      uu = rho1*(-rho2)/(d**2*(2d0*sk*g)**3)
      vv = (2d0*rho1/xr(k,ia)+rho2)/(d*(2d0*sk*g)**2)
      ww = rho1*(rhou1-rhod1-zet*rho1)
     &     /(d**2*(2d0*sk*g)**2)
c
      call excpbec(rs,zet,tt,uu,vv,ww,ec,vcup,vcdn)
c
      exc = exc + ec
      vxcup = vxcup + vcup
      vxcdn = vxcdn + vcdn
      endif
c --- convert from h to ry
      exc = 2d0*exc
      xu = 2d0*vxcup
      xd = 2d0*vxcdn
      v(k,ia,1) = v(k,ia,1) + xu
      v(k,ia,2) = v(k,ia,2) + xd
      u = w(k,ia,1)
      w(k,ia,1) = ro(k,ia,1)*(u+exc)
c     w(k,ia,2) = ro(k,ia,1)*u - 3d0*(exc-xu)*(ro(k,ia,1)+ro(k,ia,2))
c    &   *5d-1 - 3d0*(exc-xd)*(ro(k,ia,1)-ro(k,ia,2))*5d-1
   10 continue
      end
      subroutine excpbex(rho,s,u,v,ex,vx)
c----------------------------------------------------------------------
c  pbe exchange for a spin-unpolarized electronic system
c  k burke's modification of pw91 codes, may 14, 1996
c  modified again by k. burke, june 29, 1996, with simpler fx(s)
c----------------------------------------------------------------------
c  input rho : density
c  input s:  abs(grad rho)/(2*kf*rho), where kf=(3 pi^2 rho)^(1/3)
c  input u:  (grad rho)*grad(abs(grad rho))/(rho**2 * (2*kf)**3)
c  input v: (laplacian rho)/(rho*(2*kf)**2)
c   (for u,v, see pw86(24))
c  output:  exchange energy per electron (ex) and potential (vx)
c----------------------------------------------------------------------
c references:
c [a] j.p. perdew, k. burke, and m. ernzerhof,
c     phys. rev. lett. 77, 3865 (1996).
c [b] j.p. perdew and y. wang,
c     phys. rev. b33, 8800 (1986); b40, 3399 (1989) (e).
c----------------------------------------------------------------------
c formulas:
c   	e_x[unif]=ax*rho^(4/3)  [lda]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[pbe]=e_x[unif]*fxpbe(s)
c	fxpbe(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13)
c----------------------------------------------------------------------
c
c  all input and output is in atomic units!
c
c  modifications by: e. engel
c  last revision:    may 9, 2001
cengel
      implicit real*8(a-h,o-z)
c
c*** start of declarations rewritten by spag
c
c parameter definitions
c
      real*8 thrd,thrd4,ax,um,uk,ul
      parameter (thrd=1.d0/3.d0,thrd4=4.d0/3.d0,
     &           ax=-0.738558766382022405884230032680836d0,
     &           um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c
c dummy arguments
c
      real*8 ex,rho,s,u,v,vx
c
c local variables
c
      real*8 exunif,fs,fss,fxpbe,p0,s2
c
c*** end of declarations rewritten by spag
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct lda exchange energy density
      exunif = ax*rho**thrd
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct pbe enhancement factor
      s2 = s*s
      p0 = 1.d0 + ul*s2
      fxpbe = 1d0 + uk - uk/p0
      ex = exunif*fxpbe
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c  energy done. now the potential:
c  find first and second derivatives of fx w.r.t s.
c  fs=(1/s)*d fxpbe/ ds
c  fss=d fs/ds
      fs = 2.d0*uk*ul/(p0*p0)
      fss = -4.d0*ul*s*fs/p0
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c calculate potential from [b](24)
      vx = exunif*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
      end
      subroutine excpbec(rs,zeta,t,uu,vv,ww,ec,vcup,vcdn)
cengel
c  this subroutine evaluates the correlation energy per particle and
c  spin-up and spin-dn correlation potentials within the perdew-burke-
c  ernzerhof gga. it is a slightly modified version of k. burke's
c  official pbe subroutine.
c
c  input:  rs   = wigner-seitz radius = ( 3 / (4*pi*(dup+ddn)) )**(1/3)
c          zeta = relative spin polarization = (dup-ddn)/(dup+ddn)
c          t    = abs(grad d) / ( (2*sk*g) * d )
c          uu   = (grad d)*grad(abs(grad d)) / ( (2*sk*g)**3 * d**2 )
c          vv   = (laplacian d) / ( (2*sk*g)**2 * d )
c          ww   = (grad d)*(grad zeta) / ( (2*sk*g)**2 * d )
c  where:  fk   = local fermi momentum = (3*pi**2*(dup+ddn))**(1/3)
c          sk   = local screening momentum = (4*fk/pi)**(1/2)
c
c  output: ec   = correlation energy per particle
c          vcup = spin-up correlation potential
c          vcdn = spin-dn correlation potential
c
c  all input and output is in atomic units!
c
c references:
c [a] j.p. perdew, k. burke, and m. ernzerhof,
c     phys. rev. lett. 77, 3865 (1996).
c [b] j. p. perdew, k. burke, and y. wang,
c     phys. rev. b54, 16533 (1996).
c [c] j. p. perdew and y. wang,
c     phys. rev. b45, 13244 (1992).
c
c
c  last revision:    may 9, 2001
c  written by:       k. burke, may 14, 1996.
c  modifications by: e. engel
cengel
      implicit real*8(a-h,o-z)
c
c*** start of declarations rewritten by spag
c
c parameter definitions
c
      real*8 thrd,thrdm,thrd2,sixthm,thrd4,gam,fzz,gamma,bet,delt,eta
      parameter (thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd,
     &           sixthm=thrdm/2.d0,thrd4=4.d0*thrd,
     &           gam=0.5198420997897463295344212145565d0,
     &           fzz=8.d0/(9.d0*gam),
     &           gamma=0.03109069086965489503494086371273d0,
     &           bet=0.06672455060314922d0,delt=bet/gamma,eta=1.d-12)
c
c dummy arguments
c
      real*8 ec,rs,t,uu,vcdn,vcup,vv,ww,zeta
c
c local variables
c
      real*8 alfm,alfrsm,b,b2,bec,bg,comm,ecrs,eczeta,ep,eprs,eu,eurs,f,
     &       fac,fact0,fact1,fact2,fact3,fact5,fz,g,g3,g4,gz,h,hb,hbt,
     &       hrs,hrst,ht,htt,hz,hzt,pon,pref,q4,q5,q8,q9,rsthrd,rtrs,t2,
     &       t4,t6,z4
      external excgcor2
c
c*** end of declarations rewritten by spag
c
c thrd*=various multiples of 1/3
c numbers for use in lsd energy spin-interpolation formula, [c](9).
c      gam= 2^(4/3)-2
c      fzz=f''(0)= 8/(9*gam)
c numbers for construction of pbe
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at
c          |zeta|=1.
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find lsd energy contributions, using [c](10) and table i[c].
c eu=unpolarized lsd correlation energy
c eurs=deu/drs
c ep=fully polarized lsd correlation energy
c eprs=dep/drs
c alfm=-spin stiffness, [c](3).
c alfrsm=-dalpha/drs
c f=spin-scaling factor from [c](9).
c construct ec, using [c](8)
      if ( rs.lt.3.d5 ) then
         rtrs = sqrt(rs)
         call excgcor2(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0,
     &                 0.49294d0,rtrs,eu,eurs)
         call excgcor2(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,
     &                 3.3662d0,0.62517d0,rtrs,ep,eprs)
         call excgcor2(0.0168869d0,0.11125d0,10.357d0,3.6231d0,
     &                 0.88026d0,0.49671d0,rtrs,alfm,alfrsm)
         z4 = zeta**4
         f = ((1.d0+zeta)**thrd4+(1.d0-zeta)**thrd4-2.d0)/gam
         ec = eu*(1.d0-f*z4) + ep*f*z4 - alfm*f*(1.d0-z4)/fzz
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c lsd potential from [c](a1)
c ecrs = dec/drs [c](a2)
c eczeta=dec/dzeta [c](a3)
c fz = df/dzeta [c](a4)
         ecrs = eurs*(1.d0-f*z4) + eprs*f*z4 - alfrsm*f*(1.d0-z4)/fzz
         fz = thrd4*((1.d0+zeta)**thrd-(1.d0-zeta)**thrd)/gam
         eczeta = 4.d0*(zeta**3)*f*(ep-eu+alfm/fzz)
     &            + fz*(z4*ep-z4*eu-(1.d0-z4)*alfm/fzz)
         comm = ec - rs*ecrs/3.d0 - zeta*eczeta
         vcup = comm + eczeta
         vcdn = comm - eczeta
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c pbe correlation energy
c g=phi(zeta), given after [a](3)
c delt=bet/gamma
c b=a of [a](8)
         g = ((1.d0+zeta)**thrd2+(1.d0-zeta)**thrd2)/2.d0
         g3 = g**3
         pon = -ec/(g3*gamma)
         b = delt/(exp(pon)-1.d0)
         b2 = b*b
         t2 = t*t
         t4 = t2*t2
         q4 = 1.d0 + b*t2
         q5 = 1.d0 + b*t2 + b2*t4
         h = g3*(bet/delt)*log(1.d0+delt*q4*t2/q5)
         ec = ec + h
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c energy done. now the potential, using appendix e of [b].
         g4 = g3*g
         t6 = t4*t2
         rsthrd = rs/3.d0
         gz = (((1.d0+zeta)**2+eta)**sixthm-((1.d0-zeta)**2+eta)
     &        **sixthm)/3.d0
         fac = delt/b + 1.d0
         bg = -3.d0*b2*ec*fac/(bet*g4)
         bec = b2*fac/(bet*g3)
         q8 = q5*q5 + delt*q4*q5*t2
         q9 = 1.d0 + 2.d0*b*t2
         hb = -bet*g3*b*t6*(2.d0+b*t2)/q8
         hrs = -rsthrd*hb*bec*ecrs
         fact0 = 2.d0*delt - 6.d0*b
         fact1 = q5*q9 + q4*q9*q9
         hbt = 2.d0*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
         hrst = rsthrd*t2*hbt*bec*ecrs
         hz = 3.d0*gz*h/g + hb*(bg*gz+bec*eczeta)
         ht = 2.d0*bet*g3*q9/q8
         hzt = 3.d0*gz*ht/g + hbt*(bg*gz+bec*eczeta)
         fact2 = q4*q5 + b*t2*(q4*q9+q5)
         fact3 = 2.d0*b*q5*q9 + delt*fact2
         htt = 4.d0*bet*g3*t*(2.d0*b/q8-(q9*fact3/q8)/q8)
         comm = h + hrs + hrst + t2*ht/6.d0 + 7.d0*t2*t*htt/6.d0
         pref = hz - gz*t2*ht/g
         fact5 = gz*(2.d0*ht+t*htt)/g
         comm = comm - pref*zeta - uu*htt - vv*ht - ww*(hzt-fact5)
         vcup = vcup + comm + pref
         vcdn = vcdn + comm - pref
      else
         vcup = 0.d0
         vcdn = 0.d0
      end if
      end
      subroutine excgcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
c slimmed down version of gcor used in pw91 routines, to interpolate
c lsd correlation energy, as given by (10) of
c j. p. perdew and y. wang, phys. rev. b {\bf 45}, 13244 (1992).
c k. burke, may 11, 1996.
      implicit real*8(a-h,o-z)
c
c*** start of declarations rewritten by spag
c
c dummy arguments
c
      real*8 a,a1,b1,b2,b3,b4,gg,ggrs,rtrs
c
c local variables
c
      real*8 q0,q1,q2,q3
c
c*** end of declarations rewritten by spag
c
      q0 = -2.d0*a*(1.d0+a1*rtrs*rtrs)
      q1 = 2.d0*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
      q2 = log(1.d0+1.d0/q1)
      gg = q0*q2
      q3 = a*(b1/rtrs+2.d0*b2+rtrs*(3.d0*b3+4.d0*b4*rtrs))
      ggrs = -2.d0*a*a1*q2 - q0*q3/(q1*(1.d0+q1))
      end
