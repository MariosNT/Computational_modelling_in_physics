      Program quadrature                  !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 xl,xr
      real*8 func, simpson
      external func
      integer Nh1,Nh,ih,Nstep,ifact
c
      open(10,file='quadrat.inp',status='unknown')
      open(11,file='quadrat.out',status='unknown')
c
      read(10,*)xl,xr,Nh1,ifact,Nstep
c
      Nh=Nh1
c
!-----------------------------------------------------------------------
c
      do ih=1,Nstep
      	 write(11,'(i8,F15.10)')Nh,simpson(xl,xr,Nh,func)
         Nh=Nh*ifact
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
      real*8 FUNCTION simpson(xl,xr,Nh,fun)
c
      implicit none
c
      real*8 x,h,fun,xl,xr
      integer ih,Nh
c
      h=(xr-xl)/Nh/2.d0
c
      simpson=0.d0
      x=xl+h
c
      do ih=1,Nh-1
         if (x<2.d0)   then
            simpson=simpson + h*(fun(x+h)+4.d0*fun(x)+fun(x-h))/3.d0
         endif
	 x=x+2.d0*h
      enddo
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION func(x)
c
      implicit none
c
      real*8 x
c
      func=sqrt(4.d0-x**2)
c
      return
      end
c
