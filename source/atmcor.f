      subroutine atmcor(anclr,anc,nc)
c-----------------------------------------------------------------------
c     Returns shell number belonging core states.
c     coded by H.Akai, 1979, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ncore(9),nshell(10)
      data ncore/2, 8, 8, 12, 6, 12, 6, 26, 6/
     &    ,nshell/0, 1, 3, 5, 7, 8, 10, 11, 14, 15/
      nclr=anclr
      n=0
      do 10 i=1,9
      ii=i
      m=n+ncore(i)
      if(m .gt. nclr) go to 20
   10 n=m
      ii=ii+1
   20 anc=n
      nc=nshell(ii)
      return
      end
