      function asymbl(z)
c-----------------------------------------------------------------------
c     Given an atomic number z this function returns the character of
c     the corresponding atomic symbol. 'vc' is returned when z=0.
c     coded by a.akai, 1985, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension table(0:104)
      character*2 asymbl,table
      data table/'Vc'
     &          ,'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg'
     &          ,'Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti'
     &          ,'V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As'
     &          ,'Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc'
     &          ,'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I'
     &          ,'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu'
     &          ,'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta'
     &          ,'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi'
     &          ,'Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np'
     &          ,'Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'
     &          ,'?'/
      nz=z+1d-6
      nz=min(nz,104)
      nz=max(nz,0)
      asymbl=table(nz)
      return
      end
