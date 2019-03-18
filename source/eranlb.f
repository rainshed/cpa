      subroutine eranlb(v0,v1,v2,rms,itr,damp,dr,xr,meshr,kmerr)
c-----------------------------------------------------------------------
c     Check rms error and generate new imput potential by
c     Akai-Dederichs scheme. Suited for atomic calculation.
c     coded by H.Akai, 1986, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v0(meshr),v1(meshr),v2(meshr),dr(meshr),xr(meshr)
      data ami/9d-1/ , itac/1/
      eta=1d0-damp*ami
      beta=1d0+(1d0/eta-sqrt(1d0/eta**2-1d0))**2
      call rmserb(v1,v2,rms,dr,xr,kmerr)
      if(rms .lt. 1d-30) rms=1d-30
      call biomix(v1,v2,damp,meshr)
      if(itr .le. itac) go to 30
      call biomix(v0,v2,beta,meshr)
   30 call equarr(v1,v0,meshr)
      call equarr(v2,v1,meshr)
      write(6,1000)itr,log10(rms)
 1000 format('   itr=',i3,'      rms error =',4f7.3)
      return
      end
