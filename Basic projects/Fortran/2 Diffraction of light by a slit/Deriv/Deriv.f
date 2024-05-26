      Program derivatives                 !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 x,h,a
      real*8 func,fp3,fpp3
      external func
      integer Nh,ih
c
      open(10,file='deriv.inp',status='unknown')
      open(11,file='deriv.out',status='unknown')
c
      read(10,*)x,h,a,Nh
c
!-----------------------------------------------------------------------
c
      do ih=1,Nh
	 h=h/a
         write(11,'(4F15.9)')-log10(h),func(x),
     &	                     fp3(x,h,func),fpp3(x,h,func)
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
      real*8 FUNCTION fp3(x,h,fun)
c
      implicit none
c
      real*8 x,h,fun
c
      fp3=(fun(x+h)-fun(x-h))/2.d0/h
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION fpp3(x,h,fun)
c
      implicit none
c
      real*8 x,h,fun
c
      fpp3=(fun(x+h)+fun(x-h)-2.d0*fun(x))/h**2
c
      return
      end
c
!====================================================================
c
      real*8 FUNCTION func(x)
c
      implicit none
c
      real*8 x
c
      func=dsin(x)+dcos(x)
c
      return
      end
c
!====================================================================
