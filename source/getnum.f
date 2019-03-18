      subroutine getnum(ise,r,n)
c-----------------------------------------------------------------------
c     portable randum number generator
c     coded by H. Akai
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (a=7d0**5, scale=2d0**31-1)
      real*8 r(n)
      data ibias/2008120806/
c     se=ise .xor. ibias
      se=ieor(ise,ibias)
      do 10 i=1,n
      se=se*a
      in=se/scale
      se=se-scale*dble(in)
   10 r(i)=se/scale
      end
