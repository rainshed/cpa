      subroutine cstatc(z,rc,ro,corlvl,nc,cnf,ams,bms,xr,msr,wk,ier,ns1)
c-----------------------------------------------------------------------
c     Calculate core state with matching boundary condition.
c     Suitable for the atomic calculation.
c     coded by H.Akai, 1986, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension z(ns1),xr(ns1),corlvl(18),l(18)
     &         ,rc(ns1),ro(ns1),wk(ns1),cnf(18),npq(18)
      data
c                1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
     &       npq/ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6, 5, 6, 7/
     &        ,l/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0/
c
     &    ,istop/30/ , tol/1d-5/ , upper/-5d-2/
      ier=0
      call clrarr(ro,msr)
      call clrarr(rc,msr)
      do 40 j=1,18
      if(abs(cnf(j)) .lt. 1d-10) go to 40
      jj=l(j)+1
      node=npq(j)-jj
      e=corlvl(j)
      emax=1d10
      emin=-1d10
      kmatch=0
      do 10 itr=1,istop
      call radnrl(e,jj,wk,kmatch,g1,g2,nn,z,ams,bms,xr,msr,ns1)
      if(nn .gt. node) go to 100
      if(nn .lt. node) go to 110
      do 120 k=1,msr
  120 wk(k)=((ams+xr(k))*wk(k))**2
      s=0d0
      do 130 k=1,msr-2,2
  130 s=s+wk(k)+4d0*wk(k+1)+wk(k+2)
      w=3d0/s/bms
      dlt=w*wk(kmatch)*(g1-g2)/(ams+xr(kmatch))**2/bms
      if(abs(dlt) .lt. tol) go to 20
      if(abs(dlt) .gt. 1d-2*abs(e)) kmatch=0
      e=e+dlt
      if(e .gt. upper) e=5d-1*(e-dlt+upper)
      go to 10
  100 emax=min(e,emax)
      e=max(emax*1.25d0,(emax+emin)*5d-1)
      kmatch=0
      go to 10
  110 emin=max(e,emin)
      e=min(emin*7.5d-1,(emax+emin)*5d-1)
      kmatch=0
   10 continue
      write(6,1000)j
 1000 format('   ***err in cstatc...no convergence for j=',i3)
      stop
c
   20 corlvl(j)=e
      if(node .ne. nn)ier=j
      w=w*cnf(j)
      if(j .gt. nc) go to 50
      if(jj .eq. 1) rc(1)=rc(1)+w*ams
      ro(1)=rc(1)
      do 30 k=2,msr
      rc(k)=rc(k)+w*wk(k)
   30 ro(k)=rc(k)
      go to 40
   50 if(jj .eq. 1) ro(1)=ro(1)+w*ams
      do 60 k=2,msr
   60 ro(k)=ro(k)+w*wk(k)
   40 continue
      do 70 k=2,msr
      rr=1d0/xr(k)**2/(ams+xr(k))
      rc(k)=rc(k)*rr
   70 ro(k)=ro(k)*rr
      end
