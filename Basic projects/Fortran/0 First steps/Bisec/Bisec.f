      Program BISEC                       !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
      real*8 xl,xr,eps
      real*8 func,zero
c
      external func
c
      open(10,file='bisec.inp',status='unknown')
      open(11,file='bisec.out',status='unknown')
c
      read(10,*)xl,xr,eps
c
!-----------------------------------------------------------------------
c
      write(11,'(2F13.10)')ZERO(xl,xr,eps,func),eps
c
      close(10)
      close(11)
c
      stop
      end
c
!=======================================================================
c
      real*8 function ZERO(xl,xr,eps,fun)
c
      implicit none
      real*8 xl,xr,xm,eps
      real*8 fun
c
100   xm=(xl+xr)/2.d0
      if((xr-xl).lt.eps) goto 200
c
      if(fun(xl)*fun(xm).lt.0.)then
         xr=xm
         goto 100
      else
         xl=xm
         goto 100
      endif
c
200   ZERO=xm
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION FUNC(x)
c
      implicit none
      real*8 x
c
      func=dexp(x+1.d0)-1.d0
c
      return
      end
c
!=======================================================================
