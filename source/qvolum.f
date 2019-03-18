      function qvolum(anclr,ityp)
c----------------------------------------------------------------------
c     Given a nuclear charge, this funcion returns guess value of the
c     atomic volume (in the a.u.).
c     If ityp=1 is chosen, the experimental values of the atomic
c     volume for pure systems are used.
c     If ityp=2 is chosed the most values are still obtained by
c     experiments, but for the systems for which the band structure
c     calculations are available the values given in the book by
c     MJW are used.
c     (Ver. MJW), coded by H.Akai, 25 Aug. 1991.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension atvol(100,2)
      real*4 atvol
c
c     Followings are the experimental data (in cc/mol).
      data (atvol(k,1),k=1,100)
     & /14.40,0.000,13.02,4.891,4.388,5.260,17.30,14.00,17.10,16.70,
     &  23.79,14.00,10.00,12.07,13.96,17.41,22.79,28.50,45.61,26.19,
     &  15.06,10.64,8.365,7.231,7.357,7.094,6.689,6.593,7.114,9.165,
     &  11.81,13.64,12.96,16.43,23.50,38.90,56.07,33.93,19.88,14.02,
     &  10.83,9.387,8.635,8.178,8.292,8.879,10.27,13.00,15.73,16.30,
     &  18.21,20.46,25.74,37.30,69.19,38.08,22.54,17.03,20.82,20.59,
     &  20.33,19.95,28.98,19.94,19.26,18.99,18.75,18.46,18.13,24.87,
     &  17.77,13.45,10.80,9.551,8.860,8.441,8.524,9.094,10.22,14.09,
     &  17.22,18.27,21.33,22.53,0.000,50.50,73.00,38.80,22.56,19.79,
     &  15.03,13.16,13.11,12.06,0.000,0.000,0.000,0.000,0.000,0.000/
c
c     some of the following data are taken from mjw (in cc/mol).
      data (atvol(k,2),k=1,100)
     & /1.774,0.000,11.81,4.723,4.388,5.260,17.30,14.00,17.10,16.70,
     &  20.37,13.22,9.793,12.07,13.96,17.41,22.79,28.50,37.65,22.31,
     &  13.65,9.640,7.587,6.643,6.249,6.531,6.241,6.269,6.892,8.502,
     &  10.71,13.64,12.96,16.43,23.50,38.90,47.49,28.73,17.54,12.85,
     &  10.63,9.118,8.608,8.327,8.467,9.114,10.55,13.22,15.99,16.30,
     &  18.21,20.46,25.74,37.30,69.19,38.08,22.54,17.03,20.82,20.59,
     &  20.33,19.95,28.98,19.94,19.26,18.99,18.75,18.46,18.13,24.87,
     &  17.77,13.45,10.80,9.551,8.860,8.441,8.524,9.094,10.22,14.09,
     &  17.22,18.27,21.33,22.53,0.000,50.50,73.00,38.80,22.56,19.79,
     &  15.03,13.16,13.11,12.06,0.000,0.000,0.000,0.000,0.000,0.000/
      nclr=anclr
      if(nclr .lt. 0 .or. nclr .gt. 100) then
      call errtrp(1,'qvolum','anclr is illegal')
      else if(nclr .eq. 0) then
c     --- if the vacancy is specified, zero atomic volume is assumed.
      qvolum=0d0
      else
c     --- 11.205828 converts the unit from cc/mol to (a.u.)**3
      qvolum=dble(atvol(nclr,ityp))*11.205828d0
      endif
      end
