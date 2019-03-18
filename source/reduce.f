      subroutine reduce(sm,ch,tm,ng,mg)
c-----------------------------------------------------------------------
c     Generate weight for each energy point from sm-matrix.
c     coded by H.Akai, 1983, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension sm(ng,mg),ch(mg,ng),tm(ng,ng)
      do 10 i=1,mg
      do 10 j=1,ng
      ch(i,j)=0d0
      do 10 k=1,ng
   10 ch(i,j)=ch(i,j)+sm(k,i)*tm(k,j)
      return
      end
