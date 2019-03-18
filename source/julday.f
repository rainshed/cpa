      function julday(mm,id,iyyy)
c-----------------------------------------------------------------------
c     In this routine JULDAY returns the Julian Day Number which begins
c     at noon of the calendar date specified by month MM, day ID, and
c     year IYYY, all integer variables. Positive year singnifies A.D.;
c     negative, B.C. Remember that the year after 1 B.C. was 1 A.D.
c     See Numerical Recipes, Chapter 1, page 10.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (igreg=15+31*(10+12*1582))
      if(iyyy .eq. 0)
     &      write(*,*)' ***wrn in julday...there is no Year Zero'
      if(iyyy .le. 0) iyyy=iyyy+1
      if(mm .gt. 2) then
      jy=iyyy
      jm=mm+1
      else
      jy=iyyy-1
      jm=mm+13
      endif
      julday=int(365.25d0*dble(jy))+int(30.6001d0*dble(jm))+id+1720995
      if(id+31*(mm+12*iyyy) .ge. igreg) then
      ja=int(0.01d0*dble(jy))
      julday=julday+2-ja+int(0.25d0*dble(ja))
      endif
      end
