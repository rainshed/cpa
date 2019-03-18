      subroutine dsenum(type,anclr,anc,confu,confd,corlu,corld
     &                 ,tofu,tofd,mxl,ls)
c-----------------------------------------------------------------------
c     Make standard print out of the electronic structure calculation.
c     coded by H.Akai, 1983, Juelich
c     Adapted to the spin-orbit version, 20 April, 1997, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 tofu(mxl**2),tofd(mxl**2),corlu(18),corld(18)
     &      ,confu(18),confd(18),ul(6),dl(6)
      integer lst(18)
      character asymbl*2,lsymbl(18)*2,vsymbl(6)*1,type*(*),op(18)*1
      data lsymbl/'1s','2s','2p','3s','3p','3d','4s','4p','4d','5s'
     &           ,'5p','4f','5d','6s','6p','5f','6d','7s'/
     &    ,vsymbl/'s','p','d','f','g','h'/
     &    ,zero/1d-10/
      mmxl=mxl**2
      cndn=0d0
      cnup=0d0
      orup=0d0
      ordn=0d0
      do 40 j=1,6
      ul(j)=0d0
   40 dl(j)=0d0
c     write(*,'(1x,9f12.7)')tofu
c     write(*,'(1x,9f12.7)')tofd
      do 10 l=1,mmxl
      ll=sqrt(dble(l)-0.5d0)+1
      m=l-ll*(ll-1)-1
      orup=orup+dble(m)*tofu(l)
      ordn=ordn+dble(m)*tofd(l)
      cnup=cnup+tofu(l)
      cndn=cndn+tofd(l)
      ul(ll)=ul(ll)+tofu(l)
   10 dl(ll)=dl(ll)+tofd(l)
      smnt=cnup-cndn
      orbm=orup+ordn
c     --- if spin-orbit coupling is not included, the above
c         calculation never gives the correct orbital moment
c         even though orbm has a non zero valus; non vanishing
c         value simply due to the fact the real hermonics are
c         used instead of the spherical ones when spin-orbit
c         interaction is not considered.
      if(ls .eq. 0) then
      orup=0d0
      ordn=0d0
      orbm=0d0
      endif
      ctot=anc+cnup+cndn
      write(6,1000)type,asymbl(anclr),anclr,anc
 1000 format(///t30,'*** type-',a8,1x,a2,' (z=',f5.1,') ***'
     &      /'   core charge in the muffin-tin sphere =',f10.7)
      write(6,1100)'up  ',(ul(l),'(',vsymbl(l),')',l=1,mxl)
      write(6,1100)'down',(dl(l),'(',vsymbl(l),')',l=1,mxl)
 1100 format('   valence charge in the cell (spin ',a4
     &       ,') =  ',6(f10.5,3a1))
      write(6,1200)ctot,cnup,cndn,smnt,orbm,orup,ordn
 1200 format('   total charge=',f10.5,'   valence charge (up/down)='
     &       ,2f10.5
     &       /'   spin moment=',f10.5,'  orbital moment=',f10.5
     &       /'   orbital current (up/down)=',2f10.5/)
      if(ls .eq. 1) then
      write(*,'(/a)')'   orbital occupation (spin up  )'
      do 50 l=1,mxl
      m1=(l-1)**2+1
      m2=l**2
   50 write(*,'(5x,a,x,9f10.5)')vsymbl(l),(tofu(m),m=m1,m2)
      write(*,'(/a)')'   orbital occupation (spin down)'
      do 60 l=1,mxl
      m1=(l-1)**2+1
      m2=l**2
   60 write(*,'(5x,a,x,9f10.5)')vsymbl(l),(tofd(m),m=m1,m2)
      write(*,*)
      endif
      nc=0
      do 20 j=1,18
      if(abs(confu(j)) .gt. zero) then
      nc=nc+1
      lst(nc)=j
      if(confu(j) .gt. zero) then
      op(nc)=' '
      else
      op(nc)='*'
      endif
      endif
   20 continue
      if(nc .gt. 0) then
      write(6,1300)'up  ',(corlu(lst(j)),lsymbl(lst(j)),op(j),j=1,nc)
      endif
      nc=0
      do 30 j=1,18
      if(abs(confd(j)) .gt. zero) then
      nc=nc+1
      lst(nc)=j
      if(confd(j) .gt. zero) then
      op(nc)=' '
      else
      op(nc)='*'
      endif
      endif
   30 continue
      if(nc .gt. 0) then
      write(6,1300)'down',(corld(lst(j)),lsymbl(lst(j)),op(j),j=1,nc)
      endif
 1300 format('   core level  (spin ',a4,')'
c    &       /(1x,4(2x,f10.4,' Ry(',a2,')',a1)))
     &       /(1x,3(2x,f15.7,' Ry(',a2,')',a1)))
      return
      end
