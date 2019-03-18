      subroutine readk(vkp,coa,boa,nk3,nk3x,kcrt,kcx,kblst,iread)
c-----------------------------------------------------------------------
c     This program read-in the k-points according either format (a),
c     (b) or (c).
c
c     (a) simply input all the k-point to be calculated such like
c         0.0 0.0 0.0
c         0.2 0.0 0.0
c         0.4 0.0 0.0
c         0.6 0.0 0.0
c         0.8 0.0 0.0
c         1.0 0.0 0.0
c
c
c     (b) input k-point corresponding to the Gamma, K, M, etc., with
c         the number of devision of the k-mesh between that point and
c         the previous point. For example, input data such like
c         0.0 0.0 0.0   0
c         1.0 0.0 0.0  50
c
c         means that the line connecting (0, 0, 0) and (1, 0, 0)
c         is divided into 50 pieces, and on each grid points
c         including both endpoints the Bloch spectrum function
c         will be calculated.
c
c     (c) input k-point corresponding to the Gamma, K, M, etc., with
c         the total number of k-point to be used as the first data.
c         For example, input data such like
c         101
c         0.0 0.0 0.0
c         1.0 0.0 0.0
c         1.0 0.5 0.0
c
c         means that the line connecting (0, 0, 0), (1, 0, 0), and
c         (1, 0.5, 0) is divided into more or less equidistant 100
c          pieces, and on each grid points including both endpoints
c          the Bloch spectrum function will be calculated.
c
c     In all formats, data will be read-in untill any data that
c     do not follow the above format be met. Which format should be
c     adopted is controled by a parameter iread. If iread=1,
c     format (a), iread=2 for (b), and iread=3 for (c)
c     All the data shoulbe be given in the unit of 2pi*(1/a, 1/b, 1/c),
c     i.e. actual k vector is 2pi*(kx/a, ky/b, kz/c), where a, b, c
c     are the lattice constant along x, y, and z direction.
c
c     Coded by H. Akai, 4 June 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(kbzx=25)
      real*8 vkp(3,nk3x),bzke(3,kbzx),dbz(3),s(kbzx)
      character token*80
      integer kcrt(kcx)
      logical alphac
      call clrari(kcrt,kcx)
      nk3=0
c
c     --- read-in k-point according to format (a)
      if(iread .eq. 1) then
      kblst=1
      do 10 k=1,99999
      call xtoken(token,*20)
      if(alphac(token(1:1))) go to 20
      nk3=nk3+1
      if(k .gt. nk3x) call errtrp(1,'readk','nk3 too large')
      read(token,*)vkp(1,k)
      vkp(1,k)=vkp(1,k)+1d-20
      do 30 i=2,3
      token=' '
      call xtoken(token,*40)
   40 if(alphac(token(1:1)))
     &   call errtrp(1,'readk','illegal k-value read')
   30 read(token,*)vkp(i,k)
      vkp(2,k)=vkp(2,k)/boa+1d-10
   10 vkp(3,k)=vkp(3,k)/coa+1d-10
   20 call resume
      return
c
c     --- read-in k-point according to format (b)
      else if(iread .eq. 2) then
      kcrt(1)=1
      do 50 kb=1,kbzx
      call xtoken(token,*60)
      if(alphac(token(1:1))) go to 60
      if(kb .gt. kcx) call errtrp(1,'readk','kb too large')
      read(token,*)bzke(1,kb)
      do 70 i=2,3
      call xtoken(token,*80)
   80 if(alphac(token(1:1)))
     &   call errtrp(1,'readk','illegal k-value read')
   70 read(token,*)bzke(i,kb)
      bzke(2,kb)=bzke(2,kb)/boa
      bzke(3,kb)=bzke(3,kb)/coa
      call xtoken(token,*90)
   90 if(alphac(token(1:1)))
     &   call errtrp(1,'readk','illegal division number read')
      read(token,'(i5)')kdiv
      kdiv=max(0,kdiv)
      if(kb .eq. 1 .or. kdiv .le. 1) then
      nk3=nk3+1
      if(nk3 .gt. nk3x) call errtrp(1,'readk','nk3 too large')
      kend=nk3
      kcrt(kb)=kend
      kblst=kb
      do 100 i=1,3
  100 vkp(i,kend)=bzke(i,kb)+1d-20
      else
      do 120 i=1,3
  120 dbz(i)=(bzke(i,kb)-bzke(i,kb-1))/dble(kdiv)
      nk3=nk3+kdiv
      if(nk3 .gt. nk3x) call errtrp(1,'readk','nk3 too large')
      do 130 k=kend+1,nk3
      kk=k-kend
      do 130 i=1,3
  130 vkp(i,k)=vkp(i,kend)+dbz(i)*dble(kk)+1d-20
      kend=nk3
      kcrt(kb)=kend
      kblst=kb
      endif
   50 continue
   60 call resume
c     do 141 k=1,nk3
c 141 write(*,*)(vkp(1,k),i=1,3)
c     write(*,*)(kcrt(k),k=1,kblst)
      return
c
c     --- read-in k-point according to format (c)
      else if(iread .eq. 3) then
      call xtoken(token,*150)
      if(alphac(token(1:1))) go to 150
      read(token,*)nk3
      if(nk3 .gt. nk3x) call errtrp(1,'readk','nk3 too large')
      do 160 kb=1,kbzx
      call xtoken(token,*170)
      if(alphac(token(1:1))) then
      call resume
      go to 170
      endif
      if(kb .gt. kcx) call errtrp(1,'readk','kb too large')
      read(token,*)bzke(1,kb)
      do 180 i=2,3
      call xtoken(token,*190)
  190 if(alphac(token(1:1)))
     &   call errtrp(1,'readk','illegal k-value read')
  180 read(token,*)bzke(i,kb)
      bzke(2,kb)=bzke(2,kb)/boa
      bzke(3,kb)=bzke(3,kb)/coa
      kblst=kb
  160 if(kblst .gt. kbzx)
     &    call errtrp(1,'readk','too many terminal points')
  170 s(1)=0d0
      do 200 k=2,kblst
  200 s(k)=s(k-1)+sqrt((bzke(1,k)-bzke(1,k-1))**2+((bzke(2,k)-
     &           bzke(2,k-1)))**2+((bzke(3,k)-bzke(3,k-1)))**2)
c     write(*,'(/(i3,3f12.7,2x,f12.7))')
c    &      (k,(bzke(i,k),i=1,3),s(k),k=1,kblst)
      ds=s(kblst)/dble(nk3-1)
      kcrt(1)=1
      kdold=0
      do 210 i=1,3
  210 vkp(i,kcrt(1))=bzke(i,1)+1d-10
      do 220 k=2,kblst
      kdnew=s(k)/ds+5d-1
      kdiv=kdnew-kdold
      kdold=kdnew
      kcrt(k)=kcrt(k-1)+kdiv
      do 230 i=1,3
  230 dbz(i)=(bzke(i,k)-bzke(i,k-1))/dble(kdiv)
      do 220 kk=1,kdiv
      do 220 i=1,3
  220 vkp(i,kcrt(k-1)+kk)=bzke(i,k-1)+dbz(i)*dble(kk)+1d-20
  150 call resume
      return
c
c     --- illegal iread
      else
      call errtrp(1,'readk','illegal value of iread')
      endif
      end
