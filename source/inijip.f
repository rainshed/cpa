      subroutine inijip(ncmp,ntyp,ncmax)
c---------------------------------------------------------------------
c     This program gives mapping (i,j)->jip, where i indicates the
c     the type and j indicates the component.
c     Cade by H. Akai, 2 Dec. 1997, Osaka
c     Small modification by H. Akai, 31 July, 2015, Tokyo
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer ncmp(ntyp)
      integer,allocatable::ijtab(:,:),ityp(:),icmp(:)
      save ijtab,ityp,icmp
      ntc=ncmax*ntyp
      allocate(ijtab(ntyp,ncmax),ityp(ntc),icmp(ntc))
      ij0=0
      do 10 i0=1,ntyp
      do 10 j0=1,ncmp(i0)
      ij0=ij0+1
      ityp(ij0)=i0
      icmp(ij0)=j0
   10 ijtab(i0,j0)=ij0
      return
c---------------------------------------------------------------------
c     call jip(i,j,ij) gives mapping (type-i, component-j)->ij
c---------------------------------------------------------------------
      entry jip(i,j,ij)
      ij=ijtab(i,j)
      return
c---------------------------------------------------------------------
c     call jipinv(ij,i,j) gives mapping ij->(type-i, component-j)
c---------------------------------------------------------------------
      entry jipinv(ij,i,j)
      i=ityp(ij)
      j=icmp(ij)
      return
c---------------------------------------------------------------------
c     call endjip deallocate memory regions.
c---------------------------------------------------------------------
      entry endjip
      deallocate(ijtab,ityp,icmp)
      end
