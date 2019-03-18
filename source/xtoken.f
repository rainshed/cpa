      subroutine xtoken(token,*)
c----------------------------------------------------------------------
c     This program gets character string by calling 'console' program
c     and returns the token one by one. each token in the string
c     must be separated by spaces or a single comma. Continued
c     commas are understood as a blank token. Return 1 occurs when
c     end of file is detected.
c     coded by H.Akai, Aug. 1991, Juelich
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(length=1024)
      character buffer*(length),token*(*)
      save buffer,istr,idlm
      data buffer/' '/,istr/0/,idlm/length/
      token=' '
      if(idlm .ge. length .or. buffer(idlm+1:length) .eq. ' ') then
      call cnsole(5,buffer,*50)
      idlm=0
      endif
      do 10 i=idlm+1,length
      istr=i
      if(buffer(i:i) .ne. ' ') go to 20
   10 continue
   20 do 30 i=istr,length
      idlm=i
      if(buffer(i:i) .eq. ',') go to 40
      idlm=i-1
      if(i .gt. istr) then
      if(buffer(i-1:i-1) .eq. ' '
     &   .and. buffer(i:i) .ne. ' ') go to 40
      endif
   30 continue
      idlm=length+1
   40 token=buffer(istr:idlm-1)
      return
   50 istr=length+1
      return 1
c
c     --- resume the previous token.
      entry resume
      idlm=istr-1
      end
