      subroutine corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
c-----------------------------------------------------------------------
c     uniform-gas correlation of perdew and wang 1991
c     input: seitz radius (rs), relative spin polarization (zet)
c     output: correlation energy per electron (ec), up- and down-spin
c     potentials (vcup,vcdn), derivatives of ec wrt rs (ecrs)
c     & zet (eczet)
c     output: correlation contribution (alfc) to the spin stiffness
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      data gam,fzz/0.5198421d0,1.709921d0/
      thrd=1d0/3d0
      thrd4=4d0/3d0
      f=((1d0+zet)**thrd4+(1d0-zet)**thrd4-2d0)/gam
      call gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0
     &    ,0.49294d0,1d0,rs,eu,eurs)
      call gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0
     &    ,0.62517d0,1d0,rs,ep,eprs)
      call gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0
     &    ,0.49671d0,1d0,rs,alfm,alfrsm)
c     ---alfm is minus the spin stiffness alfc
      alfc=-alfm
      z4=zet**4
      ec=eu*(1d0-f*z4)+ep*f*z4-alfm*f*(1d0-z4)/fzz
c     ---energy done. now the potential:
      ecrs=eurs*(1d0-f*z4)+eprs*f*z4-alfrsm*f*(1d0-z4)/fzz
      fz=thrd4*((1d0+zet)**thrd-(1d0-zet)**thrd)/gam
      eczet=4d0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu
     &        -(1d0-z4)*alfm/fzz)
      comm=ec-rs*ecrs/3d0-zet*eczet
      vcup=comm+eczet
      vcdn=comm-eczet
      end
