      subroutine drvmsh(ef,er,ed,eb,e,mse,ids)
c---------------------------------------------------------------------
c      This program drives either 'cemesr' or 'cemesh' depending on
c      if the DOS display is required(ids=1) or not(ids=0).
c     coded by H.Akai, May 15, 1993, Osaka
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 e(mse)
      if(ids .eq. 1 .or. ids .eq. 3 .or. ids .eq. 4) then
      call cemesr(ef,er,ed,eb,e,mse)
      else
      call cemesh(ef,er,ed,eb,e,mse)
      endif
      end
