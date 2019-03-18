      subroutine guesse(anclr,corlvl)
c-----------------------------------------------------------------------
c     Gives starting enegy level for atomic calculations.
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension corlvl(18),a(15),b(15)
      data a/2.12739d0, 2.42754d0, 2.58383d0, 2.90646d0, 3.31490d0,
     &       4.08393d0, 3.47455d0, 4.02567d0, 5.14513d0, 4.12181d0,
     &       5.40802d0, 9.34949d0, 1.15378d1, 9.64960d0, 8.84248d0/
     &    ,b/1.38386d0, 4.82055d0, 5.65256d0, 1.24410d1, 1.59203d1,
     &       2.26399d1, 2.64046d1, 3.18513d1, 4.17587d1, 4.58410d1,
     &       5.67580d1, 6.45300d1, 8.48715d1, 7.83582d1, 8.85216d1/
      do 10 i=1,15
      corlvl(i)=-(anclr/b(i))**a(i)
   10 if(corlvl(i) .gt. 0d0) corlvl(i)=0d0
      do 20 i=16,18
   20 corlvl(i)=-1d0
      return
      end
