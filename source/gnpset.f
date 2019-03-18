      subroutine gnpset(cornr,nc,emx,gpt,ngpt,np,nr)
c-----------------------------------------------------------------------
c     This program arranges g vectors such that those belonging to the
c     preferred set comes first in 'gpt'. 'np' returns the number of
c     these preferred vectors.
c     coded by M. and H.Akai, 1980, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 gpt(3,ngpt),cornr(3,nc),v(3)
      np=0
      do 10 i=1,ngpt
      do 20 j=1,nc
      do 30 l=1,3
   30 v(l)=gpt(l,i)-cornr(l,j)
      if(v(1)**2+v(2)**2+v(3)**2 .lt. emx) then
      np=np+1
      do 40 l=1,3
      a=gpt(l,np)
      gpt(l,np)=gpt(l,i)
   40 gpt(l,i)=a
      go to 10
      endif
   20 continue
   10 continue
      nr=ngpt-np
      call gpsort(gpt,np)
      call gpsort(gpt(1,np+1),nr)
      end
