      function rnmrdt(a)
c----------------------------------------------------------------------
c     Given a character string, this program returns the first
c     numeric data contained in the string. If the data does not
c     contain any feasible numerical data, the program returns zero.
c     Coded by H. Akai, 8 June 2016, Tokyo
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*1 a*(*),b*256,c*256,x,s,t,u
      logical e,d,p,n,np,nd,ne,c1,c2,c3
      e(x)=x .eq. 'e' .or. x .eq. 'd' .or. x .eq. 'E' .or. x .eq. 'D'
      n(x)=x .ge. '0' .and. x .le. '9'
      p(x)=x .eq. '+' .or. x .eq. '-'
      d(x)=x .eq. '.'
      np(x)=n(x) .or. p(x)
      nd(x)=n(x) .or. d(x)
      ne(x)=n(x) .or. e(x)
      m=len(a)
      b=' '
      b(2:m+1)=a
      c=b
      rnmrdt=1d30
      do 10 i=2,m+1 
      s=b(i-1:i-1)
      t=b(i:i)
      u=b(i+1:i+1)
      c1=nd(s) .and. e(t) .and. np(u)
      c2=(np(s) .or. ne(u)) .and. d(t)
      c3=(e(s) .or. nd(u)) .and. p(t)
   10 if(.not. (n(t) .or. c1 .or. c2 .or. c3)) b(i:i)=' '
      read(b,*,err=20,end=20)rnmrdt
      do 30 i=1,m+1
   30 if(b(i:i) .eq. c(i:i)) c(i:i)=' '
      a=c(2:m+1)
   20 return
      end 
