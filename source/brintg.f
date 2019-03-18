      subroutine brintg(e,emx,redc,j,v,dr,xr,meshr)
c-----------------------------------------------------------------------
c     this program calculate breit integral needed to obtain
c     hyperfine filed.
c     coded by h.akai 1983, juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (c1=274.0720442d0, c2=1d0/c1**2, c3=24d0/9d0
     &          ,c4=19d0/9d0, c5=5d0/9d0, c6=1d0/9d0 )
      dimension v(meshr),dr(meshr),xr(meshr),r(3),h(10),dp(3),dq(3)
     &         ,st(3),ext(3)
      data c/c1/ , itrmx/5/
      crit=min(emx,0d0)
      l=j-1
      bl=dble(l*j)
      es=v(meshr)+e
      z=-xr(1)*v(1)*5d-1
c
c     ---- expansion of type p=x**ss*(1+ap*x+...) ----
      if(z .gt. 1d0) then
      gm=2d0*z*c2
      ss=sqrt(bl+1d0-2d0*z*gm)
      ap=(ss-1d0-2d0*z*gm)/(2d0*ss+1)/gm
      else
      gm=0d0
      ss=dble(l+1)
      ap=-z/ss
      endif
      do 610 k=1,3
      p0=xr(k)**ss
      p=p0*(1d0+ap*xr(k))
      q=p0*(ss-1d0+ap*ss*xr(k))/(xr(k)+gm)
      ek2=(v(k)-es)*xr(k)
      rm=xr(k)-ek2*c2
      ek2=bl/rm+ek2
      delt=dr(k)/xr(k)
      dp(k)=(p+rm*q)*delt/c3
      dq(k)=(ek2*p-q)*delt/c3
      r(k)=(p**2*(1d0+bl/(rm*c)**2)+(q/c)**2)/xr(k)**2
  610 st(k)=p*q*dr(k)/xr(k)**2
      s1=5d-1*st(1)+2d0*st(2)+st(3)
c
c     ---- adams-multon for sra case ----
c              ( forward )
      do 330 k=4,meshr-1
      ek1=ek2
      ek2=(v(k)-es)*xr(k)
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
      if(e .gt. crit) go to 500
      if(k .ge. meshr-4) go to 40
      if(ek2-ek1 .gt. 0d0 .and. ek2 .gt. 0d0) go to 40
  500 ps=p*q*dr(k)/xr(k)**2
      s1=s1+ps
  330 if(mod(k,2) .eq. 0) s1=s1+ps
c
      s1=s1-5d-1*ps
      s=2d0*s1/3d0
      go to 510
   40 kmatch=k
      pmatch=p
      p1=0d0
      q1=1d-20
      do 150 k=kmatch+3,meshr-1,4
      kback=k+1-mod(k,2)
      if(es-v(k) .lt. 7d-1*es-5d2) go to 260
  150 continue
      kback=meshr-1
      if(e .lt. -5d2) go to 260
      erel=e+(e/c)**2
      if(erel .gt. 0d0) go to 70
      er=xr(kback)*sqrt(-erel)
      h(1)=exp(-er)/er
      h(2)=-h(1)
      do 50 i=1,j
   50 h(i+2)=h(i)-dble(2*i-1)*h(i+1)/er
      p1=xr(kback)*h(j+1)
      q1=er*(h(j)-dble(j)*h(j+1)/er)-dble(l)*h(j+1)
      go to 260
   70 er=xr(kback)*sqrt(erel)
      h(1)=-sin(er)/er
      h(2)=-cos(er)/er
      do 80 i=1,j
   80 h(i+2)=-h(i)-dble(2*i-1)*h(i+1)/er
      p1=xr(kback)*h(j+1)
      q1=er*(-h(j)-dble(j)*h(j+1)/er)-dble(l)*h(j+1)
  260 ek2=(v(kback)-es)*xr(kback)
      rm=xr(kback)-ek2*c2
      ek2=bl/rm+ek2
c
c     ---- self-starting formula ----
c              ( backward )
      delt=-dr(kback)/xr(kback)
      dp1=(p1+rm*q1)*delt
      dq1=(ek2*p1-q1)*delt
      p2=p1+dp1
      q2=q1+dq1
      ek2=(v(kback-1)-es)*xr(kback-1)
      rm=xr(kback-1)-ek2*c2
      ek2=bl/rm+ek2
      delt=-dr(kback-1)/xr(kback-1)
      dp2=(p2+rm*q2)*delt
      dq2=(ek2*p2-q2)*delt
      p2=p1+5d-1*(dp1+dp2)
      q2=q1+5d-1*(dq1+dq2)
      p3=p1+2d0*dp2
      q3=q1+2d0*dq2
      do 440 itr=1,itrmx
      ek2=(v(kback-1)-es)*xr(kback-1)
      rm=xr(kback-1)-ek2*c2
      ek2=bl/rm+ek2
      delt=-dr(kback-1)/xr(kback-1)
      dp2=(p2+rm*q2)*delt
      dq2=(ek2*p2-q2)*delt
      ek2=(v(kback-2)-es)*xr(kback-2)
      rm=xr(kback-2)-ek2*c2
      ek2=bl/rm+ek2
      delt=-dr(kback-2)/xr(kback-2)
      dp3=(p3+rm*q3)*delt
      dq3=(ek2*p3-q3)*delt
      p2=p1+(5d0*dp1+8d0*dp2-dp3)/12d0
      q2=q1+(5d0*dq1+8d0*dq2-dq3)/12d0
      p3=p1+(dp1+4d0*dp2+dp3)/3d0
  440 q3=q1+(dq1+4d0*dq2+dq3)/3d0
c
      dp(1)=dp1/c3
      dp(2)=dp2/c3
      dp(3)=dp3/c3
      dq(1)=dq1/c3
      dq(2)=dq2/c3
      dq(3)=dq3/c3
      p=p3
      q=q3
      st1=p1*q1*dr(kback)/xr(kback)**2
      st2=p2*q2*dr(kback-1)/xr(kback-1)**2
      st3=p3*q3*dr(kback-2)/xr(kback-2)**2
      s2=5d-1*st1+2d0*st2+st3
c
c     ---- adams-multon for sra ----
c             ( backward )
      do 570 kk=4,kback-kmatch+1
      k=kback+1-kk
      ek2=(v(k)-es)*xr(k)
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
      ps=p*q*dr(k)/xr(k)**2
      s2=s2+ps
  570 if(mod(k,2) .eq. 0) s2=s2+ps
c
      red=(pmatch/p)**2
      if(s2*red/s1 .lt. 1d-18) red=0d0
      s=2d0*(s1+s2*red)/3d0
  510 s=-2d0*s
      call extorg(r0,r,xr)
      do 520 k=1,2
  520 ext(k+1)=-2d0*(st(k)+st(k+1))/2d0
      ext(1)=s
      do 530 k=2,3
  530 ext(k)=ext(k-1)-ext(k)
      call extorg(ext0,ext,xr)
      redc=ext0/r0
      return
      end
