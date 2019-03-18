      subroutine gtchst(eta,vc,nkp,vectk,lmxtyp,np,ngpt,gpt,tch,prr
     &                 ,hh,e1,e2,ng,mch,rpt,atmicp,nrpt,pexf,iblk
     &                 ,natm,ndmx,inv,nd,itype,lmxblk,itblk,its)
c-----------------------------------------------------------------------
c     Calculate tchebyceff expansion of the structural green's function
c     for a complex lattice. For each independent block of the
c     matrix this program call 'tchstr' such that only the necessary
c     part of the matrix is to be calculated.
c     Coded by H.Akai, 1992, Osaka
c     Variable 'mxl' extension coded (lots work!) 
c     by H. Akai, 25 Aug. 1999, Osaka
c     Modified to adapt OpenMP version by H. Akai, 6 Jan 2016, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 tch(mch,ng,nkp,ndmx),pexf(np,nkp,ndmx)
      real*8 atmicp(3,natm),dr(3),rpt(3,nrpt),gpt(3,ngpt)
     &      ,vectk(3,nkp),prr(np,nkp),hh(np,mch,nkp),vect(3)
      real*8,allocatable::drstr(:,:),h(:)
      integer iblk(natm,natm),lmxtyp(*),itype(natm),lmxblk(ndmx)
     &       ,itblk(5,natm**2)
      logical flg
      data zero/1d-15/
      allocate(drstr(3,ndmx),h(mch))
      nd=0
c     --- Check for each atom...
      call clrari(lmxblk,ndmx)
      jmx=0
      do 10 ia=1,natm
      lmxa=lmxtyp(itype(ia))
      jmx=max(jmx,2*lmxa)
      do 10 ib=1,natm
      lmxb=lmxtyp(itype(ib))
      lmx=lmxa+lmxb
      do 20 i=1,3
   20 dr(i)=atmicp(i,ia)-atmicp(i,ib)
c     --- Check if the equivalent block has already been listed.
      do 30 j=1,nd
      if((dr(1)-drstr(1,j))**2+(dr(2)-drstr(2,j))**2
     &   +(dr(3)-drstr(3,j))**2  .lt. zero) then
      iblk(ia,ib)=j
c     --- This combination may require bigger l=l_a+l_b than previously
c         assigned.
      if(lmx .gt. lmxblk(j)) lmxblk(j)=lmx
      go to 10
      endif
   30 continue
c     --- Since no equivalent block is found, this defines a new block.
c     --- However, if the system has the inversion symmetry, it is
c         expected that g(i,j)=conjg(g(j,i)) holds. This fact may be
c         used to save memory space and a little bit cpu time.
      if(ia .lt. ib .and. inv .eq. 1) then
c     --- iblk(ib,ia) must have already been determined since ia<ib.
c         Note that the index ib forms an inner loop than ia.
      iblk(ia,ib)=-iblk(ib,ia)
      go to 10
      endif
      nd=nd+1
c     --- The number of bloks now expands up to nd, but it should be
c         smaller than or equal to ndmx.
      if(nd .gt. ndmx)call errtrp(1,'gtchst','table overflows')
      iblk(ia,ib)=nd
      lmxblk(nd)=lmx
      do 40 i=1,3
   40 drstr(i,nd)=dr(i)
   10 continue
c
c     --- Then, 'iblk' is revised so that it can discriminates blocks
c         that have different angular momenta. In this table two
c         different pairs (ia,ib) and (ic,id) are classified as
c         'the same block' as long as both l_a<=l_c and l_b<=l_d,
c         or vice versa, hold.
c     --- The table is empty before starting.
      its=0
c     --- Again for each atom...
      do 50 ia=1,natm
      lmxa=lmxtyp(itype(ia))
      do 50 ib=1,natm
      lmxb=lmxtyp(itype(ib))
c     --- Reset the flag before starting.
      flg=.false.
c     --- Compare angular momenta with a table which is just being
c         constructed step by step. The table is now filled up to
c         the its-th row.
      do 60 ick=1,its
      if(iblk(ia,ib) .eq. itblk(1,ick)) then
c     --- Positional relation between a and b atoms coincides with
c         a tabulated one.
      if(lmxa .le. itblk(2,ick) .and. lmxb .le. itblk(3,ick)) go to 50
c     --- Since both l_a and l_b are smaller than the tabulated ones,
c         the remaining steps are skipped.
      if(lmxa .ge. itblk(2,ick) .and. lmxb .ge. itblk(3,ick)) then
c     --- Since both l_a and l_b are larger than the tabulated ones,
c         the table must be revised in the following way.
      if(flg) then
c     --- The table is already revised. This entry points a block that
c         is ompletly included in a bigger block, meaning that the
c         entry does not represent a unique block any more. It should
c         be discarded, and this is done by giving 0 to the first
c         column of the entry.
      itblk(1,ick)=0
      else
c     --- In this case, the entry should be revised. This is done simply
c         by increasing l_a and l_b and replace the atoms that represent
c         this block. To notify that the table is already revised, the
c         flag is set.
      flg=.true.
      itblk(2,ick)=lmxa
      itblk(3,ick)=lmxb
      itblk(4,ick)=ia
      itblk(5,ick)=ib
      endif
      endif
      endif
   60 continue
      if(.not. flg) then
c     --- Since flag has not turned true, the conditions considered in
c         the above have never been met. This means, the new entry
c         should be created in the table.
      its=its+1
c     --- The first column of the valid entry is the block number that
c         designates the equivalent block of the structural Green's
c         function.
      itblk(1,its)=iblk(ia,ib)
      itblk(2,its)=lmxa
      itblk(3,its)=lmxb
      itblk(4,its)=ia
      itblk(5,its)=ib
      endif
   50 continue
c
c     --- Since the table now may contain some void entry, which is
c         denoted by 0 of the first column of the entry, it is
c         compressed.
      ip=0
      do 70 i=1,its
      if(itblk(1,i) .gt. 0) then
      ip=ip+1
      do 80 j=1,5
   80 itblk(j,ip)=itblk(j,i)
      endif
   70 continue
c     --- The number of entries is given by its.
      its=ip
c     write(*,'(1x,''itblk='',5i3)')((itblk(i,j),i=1,5),j=1,its)
c
c     --- Now, 'iblk' are revised so as consistent with 'itblk'.
      do 90 ia=1,natm
      lmxa=lmxtyp(itype(ia))
      do 90 ib=1,natm
      lmxb=lmxtyp(itype(ib))
      flg=.false.
      do 100 ick=1,its
      if(iblk(ia,ib) .eq. itblk(1,ick)) then
      if(lmxa .le. itblk(2,ick) .and. lmxb .le. itblk(3,ick)) then
      iblk(ia,ib)=ick
c     write(*,'(1x,''iblk='',3i3,'' lmx='',2i3)')ia,ib,ick,lmxa,lmxb
      go to 90
      endif
      endif
  100 continue
c     --- The end of the above loop means the pair (ia,ib) does
c         not fit any tabulated entries. This obviously implies
c         the failure of the algorithm.
      call errtrp(1,'gtchst','constructing tables fails')
   90 continue
c     --- store (k+g)**2 and  yl for preferred set
      do 110 kp=1,nkp
      do 110 i=1,np
      do 120 j=1,3
  120 vect(j)=vectk(j,kp)+gpt(j,i)
      prr(i,kp)=vect(1)**2+vect(2)**2+vect(3)**2
      call realh(vect,h,jmx)
      do 110 ml=1,(jmx+1)**2
  110 hh(i,ml,kp)=h(ml)
c
c     --- finally, the expansion coefficients are calculated for
c         all the blocks that have different structural Green's
c         functions.
!$omp parallel do default(shared)
      do 130 i=1,nd
  130 call tchstr(eta,vc,nkp,vectk,lmxblk(i),np,ngpt,gpt,tch(1,1,1,i)
     &           ,e1,e2,ng,mch,rpt,drstr(1,i),nrpt,pexf(1,1,i))
!$omp end parallel do
      end
