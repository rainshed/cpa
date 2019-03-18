      subroutine bzmesh(nf,nk,g,r,protat,isymop,iwk,ls,fbz)
c-----------------------------------------------------------------------
c     This program generates BZ mesh.
c     Adapted to any Bravai lattices. The 1st BZ is now taken as
c     a parallelepiped. Irreducible parts of BZ is chosen by applying
c     rotations compatible with the crystal point symmetries. In the
c     same time the weithht for each k-point is also determined.
c     The inversion transposes the Green fucntion with trivial factors
c     (-1)**(l1-l2), since it simply changes R1-R2 into R2-R1.
c     It, therefore, is not necessary to consider both k and -k,
c     irrespective of its symmetry, namely, whether it has the
c     inversion symmetry or not.
c     For this reason, always a half of the BZ is considered in this
c     KKR package. The remainder of the BZ is taken into account
c     in the subroutine "bzmrot" by taking the tranpose.
c     Original version is coded by H.Akai, 1989, Osaka
c     Completely revised on 12 Jan 1997, Osaka
c
c     Aditional techniue is used to save the storage iwk.
c     Normaly iwk(2*nf,2*nf,2*nf) might be used. However, the
c     during the search of the irreducible k-point, the number of
c     k-point that are to be examined can become huge, which
c     certainly may cause the trouble. In order to avoid this
c     problem third index is assigned to each bit so as to
c     reduce the working area by factor of 31.
c
c     Again completeily revised by H.Akai, 22 Aug. 1997, Duisburg
c     Adapted to the program that take into account of both the
c     time reversal and inversion symmetries.
c     Revised by H. Akai, 16 Aug. 2009
c
c                             Time reversal symmetry
c                             No            Yes
c     inversion   No          Non           k -> -k
c     symmetry    Yes         R -> -R       k -> -k / R -> -R
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(zero=1d-20, tiny=1d-10, ibit=31, nk1x=10000)
      real*8 g(3,3),r(3,3),vkp(3,nkmx),wtkp(nkmx),protat(9,24,2)
     &      ,t(3),s(3),bc(3,3),bctyp(3)
      integer isymop(24,2),iwk(0:8*nf**3/ibit),nc(0:ibit-1)
      integer*1,allocatable::iwt(:)
      logical fbz
      save iwt,nc,dw,nfa,nfb,nfc,da,db,dc
c     --- nf=0 is for some test calculation or for the case that
c     the unit-cell contained so many atoms. Only the Gamma-point
c     is used. 'zero' is used in order to avoid ill behaviors
c     that may happen virtually during the calculation of the
c     structural Green function.
      if(nf .eq. 0) then
      nk=1
      return
      endif
c
      allocate(iwt(nk1x))
c  
c     --- if ls=0, the time reversal symmetry exists,
c         while if isymop(1,2)=1, the inversion symmetry exists.    
      if(ls .eq. 1 .and. isymop(1,2) .eq. 0) then
c     --- no way, both k and -k should be taken into account.
      invx=1
      else
c     --- -k can be elimitanted irrespective of the existence
c         of the inversion symmetry because of the time reversal symmetry
      invx=2
      endif
c     --- special treatment to include both k and -k irrespective of
c         the real crystal symmetry.
c     invx=1
c     --- special treatment to exclude -k irrespective of
c         the real crystal symmetry.
c     invx=2
c     --- check if the system is bcc, bct, or bco.
      do 120 j=1,3
      do 120 i=1,3
  120 bc(i,j)=r(i,j)+r(i,mod(j,3)+1)
      do 130 j=1,3
      bctyp(j)=0d0
      do 130 i=1,3
  130 bctyp(j)=bctyp(j)+bc(i,j)*bc(i,mod(j,3)+1)
      if(abs(bctyp(1))+abs(bctyp(2))+abs(bctyp(3)) .lt. 1d-8) then
c     if(.true.) then
c     ---the lattice must be bcc, bct, or bco
c     --- Specifically in this case, nfa=nfb=nfc should be
c         fulfilled. Otherwise, the rotated mesh does not
c         fit the unrotated mesh anymore.
      nfa=nf
      nfb=nf
      nfc=nf
      else
      ga=sqrt(g(1,1)**2+g(2,1)**2+g(3,1)**2)
      gb=sqrt(g(1,2)**2+g(2,2)**2+g(3,2)**2)
      gc=sqrt(g(1,3)**2+g(2,3)**2+g(3,3)**2)
      fct=dble(nf)/(ga*gb*gc)**(1d0/3d0)
      nfa=max(fct*ga+1d-3,1d0)
      nfb=max(fct*gb+1d-3,1d0)
      nfc=max(fct*gc+1d-3,1d0)
      endif
      if(.not. (nfa .eq. nfb .and. nfb .eq. nfc))
     &   write(*,'(3x,3(a,i3))')'nfa=',nfa,'  nfb=',nfb,'  nfc=',nfc
      nc(0)=1
      do 60 i=1,ibit-1
   60 nc(i)=nc(i-1)+nc(i-1) 
c    write(*,*)ibit,nc
      call clrari(iwk,8*nfa*nfb*nfc/ibit+1)
c     --- d is the mesh width
      da=1d0/dble(2*nfa)
      db=1d0/dble(2*nfb)
      dc=1d0/dble(2*nfc)
c     --- dw is the weight for each point
      dw=1d0/dble(nfa)/dble(nfb)/dble(nfc)/8d0
      nk=0
c     --- in the following k=d*(ia*ga+ib*gb+ic*gc) is generated
c         and tried if it is suitable as an irreducible k-point.
c
c     --- First generate full BZ from g-vectors.
      do 10 jc=1,2*nfc
      ic=jc-nfc
c     --- While jc runs starting from 0 to 2*nf-1, ic first runs
c         from 0 to nf, then from -1 to -nf+1. This is done only to
c         arrange the sequence of the k-points such that it looks
c         somewhat similar to the conventional k-point meshes.
      c=dc*dble(ic)
      xc=c*g(1,3)
      yc=c*g(2,3)
      zc=c*g(3,3)
      do 10 jb=1,2*nfb
      ib=jb-nfb
      b=db*dble(ib)
      xb=xc+b*g(1,2)
      yb=yc+b*g(2,2)
      zb=zc+b*g(3,2)
      do 10 ja=1,2*nfa
      ia=ja-nfa
c     --- if it is not yet marked, it is a new k point.
      l=2*nfa*(2*nfb*(jc-1)+jb-1)+ja-1
      iadr=l/ibit
      j=mod(l,ibit)
      if(mod(iwk(iadr)/nc(j),2) .eq. 0) then
c     write(*,*)ja,jb,jc
      a=da*dble(ia)
      t(1)=xb+a*g(1,1)
      t(2)=yb+a*g(2,1)
      t(3)=zb+a*g(3,1)
      nk=nk+1
      if(nk .gt. nk1x) call errtrp(1,'bzmesh','too many k-points')
      iwt(nk)=1
      if(fbz) go to 10
c     --- Then mark all the equivalent k-points.
      do 20 i=1,2
      do 20 iop=1,24
      if(isymop(iop,i) .eq. 1) then
      call rotatm(t,s,protat(1,iop,i),1)
      ga=s(1)*r(1,1)+s(2)*r(2,1)+s(3)*r(3,1)
      gb=s(1)*r(1,2)+s(2)*r(2,2)+s(3)*r(3,2)
      gc=s(1)*r(1,3)+s(2)*r(2,3)+s(3)*r(3,3)
      iga=int(ga/da+1000.5d0)-1000
      igb=int(gb/db+1000.5d0)-1000
      igc=int(gc/dc+1000.5d0)-1000
c     --- when ls=0, inversion is assumed irrespective of actual
c         inversion symmetry. This is taken care of in "bzmrot".
      do 30 inv=1,invx
c     --- when inv=2, the point -k will be marked irrespective
c          of the value of isymop(iop,3-i) 
      ip=(-1)**(inv-1)
      jga=mod(ip*iga+9*nfa-1,2*nfa)+1
      jgb=mod(ip*igb+9*nfb-1,2*nfb)+1
      jgc=mod(ip*igc+9*nfc-1,2*nfc)+1
      if(jga .ne. ja .or. jgb .ne. jb .or. jgc .ne. jc) then
      l=2*nfa*(2*nfb*(jgc-1)+jgb-1)+jga-1
      iadr=l/ibit
      j=mod(l,ibit)
      if(mod(iwk(iadr)/nc(j),2) .eq. 0) then
c     --- Since it is not marked, we mark it.
      iwk(iadr)=iwk(iadr)+nc(j)
c     --- This k-point is not actually used for the calculation
c         but contributes to the weight.
      iwt(nk)=iwt(nk)+1
      endif
      endif
   30 continue
      endif
   20 continue
      endif
   10 continue
c     --- Check if the sum of the weights precicely gives
c         unity. If not, something wrong must have happened.
      ichkwt=0
c     write(*,'(1x,a,i4)')' nk=',nk
      do 40 k=1,nk
c     do 110 i=1,3
c 110 vkp(i,k)=-vkp(i,k)
c     write(*,'(1x,i4,3f10.5,2x,f12.7)')k,(vkp(i,k),i=1,3),wtkp(k)
   40 ichkwt=ichkwt+iwt(k)
      chkwt=dw*dble(ichkwt)
      if(abs(chkwt-1d0) .gt. tiny) then
      write(*,'(1x,a,f12.7)')' chkwt=',chkwt
      call errtrp(1,'bzmesh','weight does not sum up to unity')
      endif
      return
      entry ftchvk(vkp,wtkp,g,nkmx,iwk,nf)
c-----------------------------------------------------------------------
      if(nkmx .eq. 1) then
      do 50 i=1,3
   50 vkp(i,1)=zero
      wtkp(1)=1d0
      return
      endif
c
      n=0
      do 140 jc=1,2*nfc
      ic=jc-nfc
c     --- While jc runs starting from 0 to 2*nf-1, ic first runs
c         from 0 to nf, then from -1 to -nf+1. This is done only to
c         arrange the sequence of the k-points such that it looks
c         somewhat similar to the conventional k-point meshes.
      c=dc*dble(ic)
      xc=c*g(1,3)
      yc=c*g(2,3)
      zc=c*g(3,3)
      do 140 jb=1,2*nfb
      ib=jb-nfb
      b=db*dble(ib)
      xb=xc+b*g(1,2)
      yb=yc+b*g(2,2)
      zb=zc+b*g(3,2)
      do 140 ja=1,2*nfa
      ia=ja-nfa
      l=2*nfa*(2*nfb*(jc-1)+jb-1)+ja-1
      iadr=l/ibit
      j=mod(l,ibit)
      if(mod(iwk(iadr)/nc(j),2) .eq. 0) then
c     --- if it is not marked, it is an irreducuble k point
      a=da*dble(ia)
      t(1)=xb+a*g(1,1)
      t(2)=yb+a*g(2,1)
      t(3)=zb+a*g(3,1)
      n=n+1
      if(n .gt. nkmx) call errtrp(1,'bzmesh','too many k-points')
      vkp(1,n)=t(1)+zero
      vkp(2,n)=t(2)+zero
      vkp(3,n)=t(3)+zero
c     --- The following 3 lines should be equivalent to the
c     above 3 lines. However, the following 3 lines gives
c     slightly different results from the above. This shoul be
c     because of the truncation errors.
c     vkp(1,n)=xb+a*g(1,1)+zero
c     vkp(2,n)=yb+a*g(2,1)+zero
c     vkp(3,n)=zb+a*g(3,1)+zero
      wtkp(n)=dw*dble(iwt(n))
      endif
  140 continue
      if(n .ne. nkmx) call errtrp(1,'bzmesg','nkmx not consistent')
c     --- shift the vectors by reciprocal vectors so that their
c         length becomes as small as possible.
      do 80 k=1,nkmx
      vv=0d0
      do 70 i=1,3
      vv=vv+vkp(i,k)**2
   70 s(i)=vkp(i,k)
      do 80 l=-1,1
      do 80 m=-1,1
      do 80 n=-1,1
      if(l**2+m**2+n**2 .ne. 0) then
      tt=0d0
      do 90 i=1,3
      t(i)=s(i)+dble(l)*g(i,1)+dble(m)*g(i,2)+dble(n)*g(i,3)
   90 tt=tt+t(i)**2
      if(tt .lt. vv) then
      vv=tt
      do 100 i=1,3
  100 vkp(i,k)=t(i)
      endif
      endif
   80 continue
      deallocate(iwt)
      end
