      subroutine cstate(v,ro,rorg,corlvl,anc,dr,xr,match
     &                 ,config,esic,meshr,wk,isr,sdftyp,ebtm,asa,mxl)
c----------------------------------------------------------------------
c     Calculate core state by matching boundary condition.
c     coded by H.Akai, 1983, Juelich
c     very minor modification for message output.
c     by H. Akai, Feb 2004.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(meshr),dr(meshr),xr(meshr),ro(meshr),wk(meshr,7)
     &         ,rorg(20),npq(18),corlvl(18),config(18),eoff(4)
      integer match(18),l(18)
      logical sic,sicon,init,asa,ifkey
      character sdftyp*12
      data small/2d-2/, amagic/1d0/
      data
c                1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
     &       npq/ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6, 5, 6, 7/
     &        ,l/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0/
c
      data istop/50/, tol/1d-8/, eb/-20d0/, sic/.true./
     &    ,eoff/ 1d3, 1d3, 1d3,-20d0/
c    &    ,eoff/ 4*-1d10/
      anclr=0d0
      do 10 k=1,meshr/8,5
   10 anclr=max(anclr,-v(k)*xr(k))
      nclr=anclr*5d-1+5d-1
      call clrarr(ro,meshr)
c     write(*,'((1x,1p6e13.6))')(v(k),k=1,meshr,10)
c     write(*,'((1x,1p6e13.6))')(-5d-1*v(k)*xr(k),k=1,60)
      anc=0d0
      tint=0d0
      esic=0d0
      do 40 j=1,18
      if(abs(config(j)) .gt. 0d0) then
c     write(*,'(a,i3,f10.3)')' j,config',j,config(j)
      jj=l(j)+1
      node=npq(j)-jj
      e=corlvl(j)
c     ei=e
      eold=e
      init=match(j) .lt. 1
      do 50 k=1,meshr
   50 wk(k,3)=v(k)
      do 60 ip=1,istop
      emax=1d10
      emin=-1d10
      if(e .gt. 1d2) e=0.1d0
c     if(abs(e) .lt. 1d-30) e=-0.1d0
      do 70 itr=1,istop
c     if(j .eq. 12)
c    &write(*,'(1x,a,i2,f12.5,2i4)')'itr,e,match,nn=',itr,e,match(j),nn
      call corada(e,jj,wk(1,1),rin,match(j),g1,g2,nn,wk(1,3)
     &           ,dr,xr,meshr,isr,asa)
c     write(*,'(1x,a,i2,f12.5,2i4,1p,2e13.6,i4)')
c    &  'itr,e,match,nn=',itr,e,match(j),nn
c    &         ,emin-eb,emax-eb,node
      if(nn .gt. node) then
      emax=min(e+eb,emax)
      e=max(emax*1.25d0,(emax+emin)*5d-1)-eb
      if(init) match(j)=0
      go to 70
      endif
      if(nn .lt. node) then
      emin=max(e+eb,emin)
      e=min(emin*0.75d0,(emax+emin)*5d-1)-eb
      if(init) match(j)=0
      go to 70
      endif
      dlt=-wk(match(j),1)*(g1-g2)
c     if(j .eq. 1)
c    & write(*,'(a,i4,1p5e20.12)')'match,wk,g1,g2,dlt,e',
c    & match(j),wk(match(j),1),g1,g2,dlt,e
      if(abs(dlt) .lt. tol) go to 80
      if(dlt .gt. 0d0) emin=max(emin,e+eb)
      if(dlt .lt. 0d0) emax=min(emax,e+eb)
      e=e+dlt
   70 continue
c     if(e .lt. 0d0) then
      if(e .lt. ebtm) then
      write(*,'(a,i3,a,i3)')
     & ' ***wrn in cstate...no convergence for nclr=',nclr,' j=',j
      write(*,'(1x,6f12.5)')(wk(k,1),k=1,meshr,10)
c     write(11,'(1x,i4/(1x,1p,6e13.6))')meshr,(xr(k),k=1,meshr)
c     write(11,'(1x,i4/(1x,1p,6e13.6))')meshr,(dr(k),k=1,meshr)
c     write(11,'(1x,i4/(1x,1p,6e13.6))')meshr,(v(k),k=1,meshr)
c     stop
      else
      e=1d3
      endif
c
   80 do 90 k=1,meshr
      wk(k,7)=v(k)-wk(k,3)
      wk(k,1)=wk(k,1)/xr(k)**2
   90 wk(k,2)=wk(k,1)-1d-10
      sicon=sic .and. e .gt. eoff(jj)
      if(.not. sicon) go to 100
c     write(*,*)anclr,j,e
      call poisna(wk(1,1),wk(1,3),0d0,dr,xr,meshr)
      do 110 k=1,meshr
  110 wk(k,5)=5d-1*wk(k,3)
      if(ifkey('vbh',sdftyp)) then
      call excvbh(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('mjw',sdftyp)) then
      call excmjw(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('vwn',sdftyp)) then
      call excvwn(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('lmm',sdftyp)) then
      call exclmm(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('pym',sdftyp)) then
      call excpym(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('pyv',sdftyp)) then
      call excpyv(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('gga91',sdftyp)) then
      call excg91(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('ev',sdftyp)) then
      call excev(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('pbe',sdftyp)) then
      call excpbe(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      else
      call errtrp(1,'cstate','sdftyp '//sdftyp//' not serviced')
      endif
      if(ip .gt. 1 .and. abs(e-eold) .lt. tol) go to 170
c     enew=5d-1*(e+eold)
      eold=e
c     e=enew
      wk(meshr,3)=v(meshr)
      do 60 k=1,meshr-1
   60 wk(k,3)=v(k)-amagic*wk(k,3)
c     write(*,'(1x,a,i3,a,i3)')'   nclr=',nclr,'   j=',j
c     call errtrp(2,'cstate','no convergence')
  170 continue
  100 corlvl(j)=e
c 100 continue
      rorg(j)=0d0
c     if(j.eq.1) write(*,'(1x,a,2f12.5)')'ei,e=',ei,e
      if(e .lt. ebtm .or. jj .gt. mxl) then
      if(sicon) then
      do 120 k=1,meshr
  120 wk(k,5)=-wk(k,7)*wk(k,1)+amagic*wk(k,5)
      esic=esic-config(j)*fintgr(wk(1,5),dr,xr,meshr)
      endif
      config(j)=abs(config(j))
c     write(*,'(1x,a)')'is included'
      anc=anc+config(j)*rin
      tint=tint+config(j)*(1d0-rin)
      do 180 k=1,meshr-1
  180 ro(k)=ro(k)+config(j)*wk(k,1)*xr(k)**2
      call extorg(rorg(j),wk,xr)
      if(abs(rorg(j)) .lt. 1d-20) rorg(j)=0d0
      rorg(j)=config(j)*rorg(j)
c     write(*,'(1x,a,3f12.6)')'e,esic,rin=',e,esic,rin
      else
      config(j)=-abs(config(j))
      endif
      if(.not. sicon .and. abs(e-ebtm) .lt. small) then
c     write(*,'(a,i3,a,2f12.5)')
c    & '   ***msg in cstate...corelevel near ebtm found for nclr='
c    &     ,nclr,' e, ebtm=',e,ebtm
      write(*,'(a,i3)')
     & '   ***msg in cstate...corelevel near ebtm found for nclr=',nclr
      endif
      endif
   40 continue
      ro(meshr)=tint
      end
