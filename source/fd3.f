      function fd3(x)
c-----------------------------------------------------------------------
c     calculates d3 appearing in the structure Green function.
c     For negative E/eta=-x (x>0), following formula is used.
c         Sum { (-x)^n/n!/(2n-1) }=-sqrt(pi*x)*erf(sqrt(x))-exp(-x)
c     On the other hand, for positive E/eta=x (x>0)
c         Sum { (x)^n/n!/(2n-1) }=exp(x)(2*sqrt(x)*F(sqrt(x))-1)
c     is used. Here F(x) is Dawson's integral:
c         F(z)=exp(-z^2) Integral (0 through z) exp(t^2) dt.
c     Both the relations can be derived by use of the recurrence
c     formula for the incomplete gamma function and the expression for 
c     the error function with the incomplete gamma function. Also it
c     can be directly obtained by use of Eq. 7.1.5 and Eq. 7.4.35 of
c     Handbook of Mathematical Functions (Abramowitz and Stegun, Dover).
c
c     Coded by H.Akai, 9 April, 2002, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      pi=4d0*atan(1d0)
      if(x .lt. 0d0) then
      fd3=-sqrt(-pi*x)*(1d0-erfc(sqrt(-x)))-exp(x)
      return
      else
      fd3=exp(x)*(2d0*sqrt(x)*dawson(sqrt(x))-1d0)
      endif
      end
