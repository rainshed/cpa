      subroutine viomix(v1,v2,a,n,meshr)
c-----------------------------------------------------------------------
c     Mix v1 and v2 with different mixing parameter for charge
c     and spin components.
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v1(n,2),v2(n,2),a(2),b(2)
      b(1)=1d0-a(1)
      b(2)=1d0-a(2)
      do 10 k=1,meshr
      vp=b(1)*(v1(k,1)+v1(k,2))+a(1)*(v2(k,1)+v2(k,2))
      vm=b(2)*(v1(k,1)-v1(k,2))+a(2)*(v2(k,1)-v2(k,2))
      v2(k,1)=5d-1*(vp+vm)
   10 v2(k,2)=5d-1*(vp-vm)
      return
      end
