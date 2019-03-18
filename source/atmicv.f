      subroutine atmicv(anclr,cnf,e,ro,dr,xr,z0,z1,z2,wk,rc,msr,ier)
c-----------------------------------------------------------------------
c              ***** program header *****
c       This program calculates a self-consistent atomic
c       potential and constructes a crystal potential
c       from it according mattheiss' prescription.
c
c       coded by H.Akai on Sept.   1980
c       revised on Oct. 23,  1982  ( at Osaka)
c       revised on Aug.  3,  1983  ( at Juelich)
c       revised on Dec. 12,  1984  ( at Juelich)
c
c       new version on Jan. 10, 1986 (at Osaka)
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension cnf(18),e(18),dr(msr),xr(msr)
     &         ,npq(18),l(18),ncnf(18),z0(msr),z1(msr),z2(msr),wk(msr)
     &         ,rc(msr),ro(msr),indx(3,19)
      character*1 lsym(4)
c
c              1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
c     ------------------------------------------------------------------
c         j     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
c     ------------------------------------------------------------------
      data npq/ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6, 5, 6, 7/
     &     , l/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0/
     &   ,ncnf/ 2, 2, 6, 2, 6,10, 2, 6,10, 2, 6,14,10, 2, 6,14,10, 2/
c
     &   ,indx/1,0,0, 2,3,0, 3,0,0, 4,5,0, 5,0,0, 7,6,0, 7,8,0, 8,0,0
     &        ,10,9,0, 10,11,0, 11,0,0, 14,12,13, 13,0,0, 14,15,0
     &        ,15,0,0, 18,16,17 ,17,0,0, 18,0,0, 0,0,0/
     &   ,lsym/'s','p','d','f'/
     &   , itrstp/30/ , tol/1d-6/ , dmp/.4d0/
      call atmcor(anclr,anc,nc)
      if(anclr .lt. 1d-10) then
c     --- if the vacancy is specified, the atomic radial mesh and the
c         atomic energy level for a hydrogen atom are assumed for
c         computational convenience.
      call mshatm(ams,bms,dr,xr,1d0,msr,msr)
      call guesse(anclr,e)
      do 100 k=1,msr
      rc(k)=0d0
  100 ro(k)=0d0
      return
      endif
      nf=0
      anv=anclr
      if(nc .lt. 1) go to 90
      do 10 i=1,nc
      if(cnf(i) .gt. 0d0) go to 50
   10 continue
      do 20 i=1,nc
      anv=anv-dble(ncnf(i))
   20 cnf(i)=ncnf(i)
      if(nc .ge. 18) go to 50
   90 do 30 i=nc+1,18
      if(cnf(i) .gt. 0d0) go to 50
   30 continue
      do 40 i=1,3
      j=indx(i,nc+1)
      if(j .eq. 0) go to 50
      if(anv .lt. 0d0) go to 50
      cnf(j)=min(anv,dble(ncnf(j)))
   40 anv=anv-dble(ncnf(j))
   50 do 80 j=1,18
   80 if(cnf(j) .gt. 0d0) nf=j
      call mshatm(ams,bms,dr,xr,anclr,msr,msr)
      call clrarr(z1,msr)
      call guessz(anclr,z1,xr,msr,msr)
      call guesse(anclr,e)
      call utimer(time,0)
      do 60 itr=1,itrstp
      itr1=itr
      call cstatc(z1,rc,ro,e,nc,cnf,ams,bms,xr,msr,wk,ier,msr)
c     write(*,*)(e(i),i=1,nf)
      call poisnb(ro,z2,anclr,ams,bms,xr,msr,msr)
      do 70 k=msr,1,-1
      z2(k)=z2(k)-fldf(ro(k))*xr(k)
      if(z2(k) .lt. 2d0) kmerr=k
   70 z2(k)=max(2d0,z2(k))
      itrd=itr
      call eranlb(z0,z1,z2,rms,itrd,dmp,dr,xr,msr,kmerr)
      if(abs(rms) .lt. tol) go to 240
   60 continue
      write(6,1500)
 1500 format('   not converge')
  240 call utimer(time,itr1)
      write(6,1300)
 1300 format (/t10,'nl',t18,'cnf',t30,'energy'/t5,35('-'))
      do 270 i=1,nf
      if(abs(cnf(i)) .lt. 1d-10) go to 270
      j=l(i)+1
      write(6,1400)npq(i),lsym(j),cnf(i),e(i)
  270 continue
      write(6,1600)
 1600 format(//)
 1400 format (t10,i1,a1,t16,f6.3,4x,f12.4)
      end
