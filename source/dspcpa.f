      subroutine dspcpa(t,m)
      complex*16 t(m,m)
      write(*,'((1x,9f12.7))')(abs(t(j,j)),j=1,m)
      end
