      Program FTABLE                      !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
      real*8 xl,xr,h,x
      real*8 myfunction
      integer Np,ip
c
      open(10,file='Ftable.inp',status='unknown')
      open(11,file='Ftable.out',status='unknown')
c
      read(10,*)xl,xr,Np
c
      h=(xr-xl)/(Np-1)
c
!-----------------------------------------------------------------------
c
      do ip=1,Np
	 x=xl+(ip-1)*h
	 write(11,'(2F15.9)')x,myfunction(x)
      enddo
c
      close(10)
      close(11)
c
      stop
      end
c
!=======================================================================
c
      real*8 FUNCTION myfunction(x)
c
      implicit none
      real*8 x
c
      myfunction=x*dsin(x)
c
      return
      end
c
!=======================================================================
