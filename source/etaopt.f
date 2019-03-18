      subroutine etaopt(emx,vc,eta,rmx,gmx)
c-----------------------------------------------------------------------
c     Given 'emx' and 'vc', this program determines the optimal 'eta'
c     parameter used in the Ewalt sum appearing in the calculation of
c     the structural Green's function. Here, emx is the maximum absolute
c     value (in the unit of (2pi/a)**2) of energy used in KKR and vc is
c     the unit cell volume (in the unit of the lattice constant a**3).
c     rmx and gmx is the maximum length of the real and reciprocal
c     lattice vectors to be used in the Ewalt sum (the unit for rmx
c     is a and that for gmx is 2pi/a.) Basic equations determining
c     rmx, gmx, and eta are the following:
c
c     (1)  (gmx**2-emx)/eta=B           (convergence in G sum)
c     (2)  eta*(pi*rmx)**2-emx/eta=B    (convergence in R sum)
c     (3)  rmx**3/vc=gmx**3*vc          (equal number for G's and R's)
c     (4)  emx/eta < alim               (exp(emx/eta) not be too large)
c
c     Here B is a number for which exp(-B) is very small.
c     In order to avoid problem that occurs because of exponetial
c     increases in all D_1, D_2, and D_3 as functions of e/eta, which
c     awfully degrades the accuracy of the structural Green's funcion
c     for large e/eta, we put a limitation emx/eta < alim.
c
c     Coded by H. Akai June 5, 1993, Osaka
c     Revised by H. Akai, 6 April 1996, Osaka
c     In the previosu version only emx=0 case is considered, but this
c     causes problems when the lattice constant is very large and
c     acrodgingly emx in the unit of (2pi/a)**2 becomes large.
c     Although the value of eta is unchanged from that of emx=0 
c     the case, the values of gmx and rmx are modified.
c     Modified by H. Akai, 7 Jan 2016, Tokyo
c     Modified by H. Akai, 26 Feb. 2017, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     parameter (zero=1d-8)
c     parameter (zero=1d-19)
      parameter (zero=1d-15)
      data alim/30d0/
      pi=4d0*atan(1d0)
      b=-log(zero)
      preta=1d0/vc**(2d0/3d0)/pi
      eta=max(preta,abs(emx/alim))
      gmx=sqrt(b*eta+emx)
      rmx=gmx/eta/pi
c     --- eta_new used in "tchstr" is different from eta defined here.
c         The relation between the two is eta_new=sqrt(eta)/2
      preta=5d-1*sqrt(preta)
      eta=5d-1*sqrt(eta)
c     eta=1.3d0*eta
c     eta=5d-1*eta
      write(*,'(/a,f8.5,a,f8.5)')'   preta=',preta,'  eta=',eta
c     write(*,'(a,f8.5,a,f8.5,a,f8.5)')
c    &             '   emx  =',emx,'   rmx  =',rmx,'  gmx=',gmx
c     stop
      end
