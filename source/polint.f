      function polint(x,xa,ya,n)
c----------------------------------------------------------------------
c     Given arrays xa and ya, each of length n and given value x,
c     this function returns a value polint. If p(x) is the polynominal
c     of degree ndg such that p(xa(i))=ya(i), i=ns,..,ns+ndg then
c     the returned value polint=p(x). ns is obtained by hunting.
c     See Numerical Recipes
c     coded by H.Akai
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (ndgmx=4, nmx=ndgmx+1)
      dimension xa(n),ya(n),c(nmx),d(nmx)
      logical ascnd
      save jlo
      data jlo/0/ , small/1d-30/
      ndg=min(ndgmx,n-1)
      ndt=ndg+1
      ascnd=xa(n) .gt. xa(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
      jlo=0
      jhi=n+1
      go to 30
      endif
      inc=1
      if(x .gt. xa(jlo) .eqv. ascnd) then
   10 jhi=jlo+inc
      if(jhi .gt. n) then
      jhi=n+1
      else if(x. gt. xa(jhi) .eqv. ascnd) then
      jlo=jhi
      inc=inc+inc
      go to 10
      endif
      else
      jhi=jlo
   20 jlo=jhi-inc
      if(jlo .lt. 1) then
      jlo=0
      else if(x .lt. xa(jlo) .eqv. ascnd) then
      jhi=jlo
      inc=inc+inc
      go to 20
      endif
      endif
   30 if(jhi-jlo .ne. 1) then
      jm=(jhi+jlo)/2
      if(x .gt. xa(jm) .eqv. ascnd) then
      jlo=jm
      else
      jhi=jm
      endif
      go to 30
      endif
      nlo=max(1,jlo-ndg/2)
      nhi=min(n,nlo+ndg)
      nlo=nhi-ndg
      if(jlo .eq. 0) then
      ns=1
      else if(jlo .eq. n) then
      ns=ndt
      else if(abs(x-xa(jlo)) .lt. abs(x-xa(jhi))) then
      ns=jlo-nlo+1
      else
      ns=jhi-nlo+1
      endif
      do 40 i=1,ndt
      ii=nlo+i-1
      c(i)=ya(ii)
   40 d(i)=ya(ii)
      polint=ya(nlo+ns-1)
      ns=ns-1
      do 60 m=1,ndg
      do 50 i=1,ndt-m
      ii=nlo+i-1
      ho=xa(ii)-x
      hp=xa(ii+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
c
c     an error can occur if two xa's are identical
      if(abs(den) .lt. small) then
      write(6,1000)
 1000 format('   ***wrn in polint...data error')
      stop
      endif
c
      den=w/den
      d(i)=hp*den
   50 c(i)=ho*den
      if(2*ns .lt. ndt-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
   60 polint=polint+dy
      return
      end
