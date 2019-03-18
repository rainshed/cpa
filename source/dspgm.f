      subroutine dspgm(gm,n1,n2,ind)
      implicit real*8 (a-h,o-z)
      complex*16 gm(n1,n1,n2)
      if(ind .eq. 1 )then
      do 10 k=n2,n2
c  10 write(*,'(1p6e13.6)')((gm(i,j,k),i=1,n1),j=1,n1)
   10 write(*,'(''gm''/(1x,1p9e10.3))')((dble(gm(i,j,k)),j=1,n1),i=1,n1)
      else if(ind. eq. 2) then
      do 20 k=n2,n2
   20 write(*,'(''gm''/(1x,1p9e10.3))')((dble(gm(j,i,k)),j=1,n1),i=1,n1)
      endif
      end
