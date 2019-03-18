      subroutine radsra(e,j,r,kmatch,g1,g2,node,z,a,b,xr,meshr,ns1)
c----------------------------------------------------------------------
c     Solves schroedinger equation for atomic states.
c     modified from corada.f by H. Akai, Nov 2011
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (c1=274.0720442d0, c2=1d0/c1**2, c3=24d0/9d0
     &          ,c4=19d0/9d0, c5=5d0/9d0, c6=1d0/9d0, msr=441)
      real*8 z(meshr),dr(msr),xr(meshr),r(meshr),h(10),dp(3),dq(3)      
      logical sra,sn1,sn2,asa
      data c/c1/ , itrmx/5/ , small/1d-2/, isr/0/, asa/.true./
      sra=isr .eq. 1
      l=j-1
      bl=dble(l*j)
      es=e
      xr(1)=1d-10
      do 10 k=1,meshr
      kb=k
      if(-z(k)*xr(k)+bl .lt. 0d0) go to 20
   10 continue
c
c     ---- expansion of type p=x**(l+1)+... ---
   20 do 610 k=1,3
      dr(k)=b*(a+xr(k))
      ek2=-z(k)-es*xr(k)
      rm=xr(k)
      if(sra) rm=rm-ek2*c2
      ek2=bl/rm+ek2
      delt=dr(k)/xr(k)
      p=xr(k)**j
      q=dble(l)*p/rm
      dp(k)=(p+rm*q)*delt/c3
      dq(k)=(ek2*p-q)*delt/c3
      r(k)=p**2
  610 if(sra) r(k)=r(k)*(1d0+bl/(rm*c)**2)+(q/c)**2
      s1=5d-1*r(1)*dr(1)+2d0*r(2)*dr(2)+r(3)*dr(3)
      node=0
      sn1=p .gt. 0d0
      if(sra) go to 300
c
c     ---- adams-multon for nrl case ----
c              ( forward )
      do 30 k=4,meshr-3
      dr(k)=b*(a+xr(k))
      ek2=-z(k)-es*xr(k)
      rm=xr(k)
      ek2=bl/rm+ek2
      dlt=c3*xr(k)/dr(k)
      gm=(1d0-dlt)/ek2
      pm=p+c4*dp(3)-c5*dp(2)+c6*dp(1)
      qm=q+c4*dq(3)-c5*dq(2)+c6*dq(1)
      q=(gm*qm-pm)/(rm+gm*(dlt+1d0))
      p=(pm+rm*q)/(dlt-1d0)
      dp(1)=dp(2)
      dp(2)=dp(3)
      dq(1)=dq(2)
      dq(2)=dq(3)
      dp(3)=p+rm*q
      dq(3)=ek2*p-q
      p=p*dlt
      q=q*dlt
      sn2=p .gt. 0d0
      if(sn1 .neqv. sn2) node=node+1
      sn1=sn2
      if(kmatch .lt. 1 .and. k .gt. kb .and. ek2 .gt. 0d0) go to 40
      if(k .eq. kmatch) go to 40
      r(k)=p**2
      ps=dr(k)*r(k)
      s1=s1+ps
   30 if(mod(k,2) .eq. 0) s1=s1+ps
      go to 310
c
c
c     ---- adams-multon for sra case ----
c              ( forward )
  300 do 330 k=4,meshr-3
      dr(k)=b*(a+xr(k))
      ek2=-z(k)-es*xr(k)
      rm=xr(k)-ek2*c2
      ek2=bl/rm+ek2
      dlt=c3*xr(k)/dr(k)
      gm=(1d0-dlt)/ek2
      pm=p+c4*dp(3)-c5*dp(2)+c6*dp(1)
      qm=q+c4*dq(3)-c5*dq(2)+c6*dq(1)
      q=(gm*qm-pm)/(rm+gm*(dlt+1d0))
      p=(pm+rm*q)/(dlt-1d0)
      dp(1)=dp(2)
      dp(2)=dp(3)
      dq(1)=dq(2)
      dq(2)=dq(3)
      dp(3)=p+rm*q
      dq(3)=ek2*p-q
      p=p*dlt
      q=q*dlt
      sn2=p .gt. 0d0
      if(sn1 .neqv. sn2) node=node+1
      sn1=sn2
      if(kmatch .lt. 1 .and. k .gt. kb .and. ek2 .gt. 0d0) go to 40
      if(k .eq. kmatch) go to 40
      r(k)=p**2*(1d0+bl/(rm*c)**2)+(q/c)**2
      ps=dr(k)*r(k)
      s1=s1+ps
  330 if(mod(k,2) .eq. 0) s1=s1+ps
c
  310 s1=s1-2d0*ps
      k=meshr-3
   40 kmatch=k
      pmatch=p
      qmatch=q
c     write(*,*)p,q,kmatch
c
      xzero=sqrt(dabs(2d3/es))
      do 150 k=kmatch+3,meshr,4
      kback=k+1-mod(k,2)
      if(xr(k) .gt. xzero) go to 130
  150 continue
      kback=meshr
c
c     --- following code gives the natural boundary condition, which is
c         suitable for a single muffin-tin potential in a free space.
c         This, however, causes some problem for positive energies.
c         In order to overcome the difficulty, I introduce somewhat
c         artificcial boundary condition for positive, or negative
c         but nearly zero, energies.
  130 ee=es
      erel=ee
      if(sra) erel=ee+(ee/c)**2
      erel=5.d-1*(erel-sqrt(erel**2+small**2))
      er=xr(kback)*sqrt(-erel)
      h(1)=exp(-er)/er
      h(2)=-h(1)
      do 50 i=1,j
   50 h(i+2)=h(i)-dble(2*i-1)*h(i+1)/er
      p1=xr(kback)*h(j+1)
      q1=er*(h(j)-dble(j)*h(j+1)/er)*ee/erel
c     ----- next statement gives normalization in the entire space---
      s3=0d0
      if(kback .ge. meshr)
     &         s3=5d-1*xr(meshr)**3*(h(j)*h(j+2)-h(j+1)**2)
c     --- normalization special
      if(asa) s3=0d0
c     p1=0d0
c     q1=0d0
      ek2=-z(kback)-es*xr(kback)
      rm1=xr(kback)
      if(sra) rm1=rm1-ek2*c2
      ek2=bl/rm1+ek2
c
c     ---- self-starting formula ----
c              ( backward )
      dr(kback)=b*(a+xr(kback))
      dr(kback-1)=b*(a+xr(kback-1))
      dr(kback-2)=b*(a+xr(kback-2))
      delt=-dr(kback)/xr(kback)
      dp1=(p1+rm1*q1)*delt
      dq1=(ek2*p1-q1)*delt
      p2=p1+dp1
      q2=q1+dq1
      ek2=-z(kback-1)-es*xr(kback-1)
      rm2=xr(kback-1)
      if(sra) rm2=rm2-ek2*c2
      ek2=bl/rm2+ek2
      delt=-dr(kback-1)/xr(kback-1)
      dp2=(p2+rm2*q2)*delt
      dq2=(ek2*p2-q2)*delt
      p2=p1+5d-1*(dp1+dp2)
      q2=q1+5d-1*(dq1+dq2)
      p3=p1+2d0*dp2
      q3=q1+2d0*dq2
      do 90 itr=1,itrmx
      ek2=-z(kback-1)-es*xr(kback-1)
      rm2=xr(kback-1)
      if(sra) rm2=rm2-ek2*c2
      ek2=bl/rm2+ek2
      delt=-dr(kback-1)/xr(kback-1)
      dp2=(p2+rm2*q2)*delt
      dq2=(ek2*p2-q2)*delt
      ek2=-z(kback-2)-es*xr(kback-2)
      rm3=xr(kback-2)
      if(sra) rm3=rm3-ek2*c2
      ek2=bl/rm3+ek2
      delt=-dr(kback-2)/xr(kback-2)
      dp3=(p3+rm3*q3)*delt
      dq3=(ek2*p3-q3)*delt
      p2=p1+(5d0*dp1+8d0*dp2-dp3)/12d0
      q2=q1+(5d0*dq1+8d0*dq2-dq3)/12d0
      p3=p1+(dp1+4d0*dp2+dp3)/3d0
   90 q3=q1+(dq1+4d0*dq2+dq3)/3d0
c
      if(sra) go to 100
      r(kback)=p1**2
      r(kback-1)=p2**2
      r(kback-2)=p3**2
      go to 110
  100 r(kback)=p1**2*(1d0+bl/(rm1*c)**2)+(q1/c)**2
      r(kback-1)=p2**2*(1d0+bl/(rm2*c)**2)+(q2/c)**2
      r(kback-2)=p3**2*(1d0+bl/(rm3*c)**2)+(q3/c)**2
  110 dp(1)=dp1/c3
      dp(2)=dp2/c3
      dp(3)=dp3/c3
      dq(1)=dq1/c3
      dq(2)=dq2/c3
      dq(3)=dq3/c3
      p=p3
      q=q3
      s2=5d-1*r(kback)*dr(kback)+2d0*r(kback-1)*dr(kback-1)
     &   +r(kback-2)*dr(kback-2)
c
c     ---- adams-multon ----
c          ( backward )
      if(sra) go to 500
c
c     ---- adams-multon for nrl ----
c             ( backward )
      do 170 kk=4,kback-kmatch+1
      k=kback+1-kk
      dr(k)=b*(a+xr(k))
      ek2=-z(k)-es*xr(k)
      rm=xr(k)
      ek2=bl/rm+ek2
      dlt=-c3*xr(k)/dr(k)
      gm=(1d0-dlt)/ek2
      pm=p+c4*dp(3)-c5*dp(2)+c6*dp(1)
      qm=q+c4*dq(3)-c5*dq(2)+c6*dq(1)
      q=(gm*qm-pm)/(rm+gm*(dlt+1d0))
      p=(pm+rm*q)/(dlt-1d0)
      dp(1)=dp(2)
      dp(2)=dp(3)
      dq(1)=dq(2)
      dq(2)=dq(3)
      dp(3)=p+rm*q
      dq(3)=ek2*p-q
      p=p*dlt
      q=q*dlt
      r(k)=p**2
      ps=dr(k)*r(k)
      s2=s2+ps
  170 if(mod(k,2) .eq. 0) s2=s2+ps
      go to 510
c
c     ---- adams-multon for sra ----
c             ( backward )
  500 do 570 kk=4,kback-kmatch+1
      k=kback+1-kk
      dr(k)=b*(a+xr(k))
      ek2=-z(k)-es*xr(k)
      rm=xr(k)-ek2*c2
      ek2=bl/rm+ek2
      dlt=-c3*xr(k)/dr(k)
      gm=(1d0-dlt)/ek2
      pm=p+c4*dp(3)-c5*dp(2)+c6*dp(1)
      qm=q+c4*dq(3)-c5*dq(2)+c6*dq(1)
      q=(gm*qm-pm)/(rm+gm*(dlt+1d0))
      p=(pm+rm*q)/(dlt-1d0)
      dp(1)=dp(2)
      dp(2)=dp(3)
      dq(1)=dq(2)
      dq(2)=dq(3)
      dp(3)=p+rm*q
      dq(3)=ek2*p-q
      p=p*dlt
      q=q*dlt
      r(k)=p**2*(1d0+bl/(rm*c)**2)+(q/c)**2
      ps=dr(k)*r(k)
      s2=s2+ps
  570 if(mod(k,2) .eq. 0) s2=s2+ps
c
  510 g1=q/p
      g2=qmatch/pmatch
      red=(pmatch/p)**2
      if(s2*red/s1 .lt. 1d-18) red=0d0
      s3=s3*red
      s=2d0*(s1+s2*red)/3d0+s3
      rin=1d0-s3/s
      red1=1d0/s
      red2=red1*red
      zero=1d-36/(abs(red2)+1d-36)
      do 120 k=1,kmatch-1
  120 r(k)=r(k)*red1
      do 140 k=kmatch,kback
      if(abs(r(k)) .lt. zero) r(k)=0d0
  140 r(k)=r(k)*red2
      do 200 k=kback+1,meshr
  200 r(k)=0d0
      do 620 k=1,kback
  620 r(k)=sqrt(r(k)/dr(k))
      g1=-g1*dr(kmatch)
      g2=-g2*dr(kmatch)
      end
