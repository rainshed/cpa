      subroutine erranc(v0,v1,v2,damp,eta,beta,cr,tend,rms,cnv,dr,xr
     &                 ,meshr,ntyp)
c-----------------------------------------------------------------------
c     Using input, output and preceding potential data, this program
c     make analysis about convergence quality and construct new trial
c     potential after Tchebyceff accelaration.
c     For detail, see H.Akai and P.P.Dederichs, J.Phys.C 18(1985)2455.
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v0(meshr,ntyp,2),v1(meshr,ntyp,2),v2(meshr,ntyp,2),c(2)
     &         ,dr(meshr,ntyp),xr(meshr,ntyp),rms(ntyp,2),tend(10)
      data enh/8d0/, s,x,x2,sx/10d0,55d0,385d0,2.8722813d0/
     &    ,zero/1d-30/
      cnv=0d0
      do 10 is=1,2
      do 10 i=1,ntyp
      call rmserr(v1(1,i,is),v2(1,i,is),rms(i,is),dr(1,i)
     &           ,xr(1,i),meshr)
      if(rms(i,is) .lt. zero) rms(i,is)=zero
   10 cnv=max(cnv,rms(i,is))
      do 20 i=1,9
   20 tend(i)=tend(i+1)
      tend(10)=log(cnv)
      y=0d0
      y2=0d0
      xy=0d0
      do 30 i=1,10
      y=y+tend(i)
      y2=y2+tend(i)**2
   30 xy=xy+dble(i)*tend(i)
      p=(s*xy-x*y)/(s*x2-x**2)
      sy=sqrt(y2/s-(y/s)**2)
      cr=p*sx/sy
      p=exp(p)
      eta=(p**2+beta)/(1d0+beta)/p
      c(1)=damp
      c(2)=min(enh*damp,5d-1)
      call viomix(v1,v2,c,ntyp*meshr,ntyp*meshr)
      c(1)=1d0+beta
      c(2)=c(1)
      call viomix(v0,v2,c,ntyp*meshr,ntyp*meshr)
      call equarr(v1,v0,2*ntyp*meshr)
      call equarr(v2,v1,2*ntyp*meshr)
      return
      end
