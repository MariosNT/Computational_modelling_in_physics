      Program waveguide                   !(P.Scharoch,R.Szymon 01.2020)
c
      implicit none
c
      real*8 h,xE,x0
      real*8 K1,K2,Kmax,deltK,epsK
      real*8 ZERO,fiE,A1,A2
c
      integer Np,ip
c
      external fiE
c
      common /param/x0,h
      common /param2/Np
c
      open(10,file='waveguide.inp',status='unknown')
      open(11,file='waveguide.out',status='unknown')
c
      read(10,*)x0,Np,Kmax
      xE=1.d0
c
      h=(xE-x0)/Np
c
      deltK=Kmax/10.d0
      epsK=Kmax/10000.d0
      K1=0.001d0*Kmax
      ip=0
c
!-----------------------------------------------------------------------
c
100   K2=K1+deltK
c
      A1=fiE(K1)
      A2=fiE(K2)
      if(A1*A2.lt.0.d0)   then
         ip=ip+1
	 write(11,'(i3,F12.6)')ip,ZERO(K1,K2,epsK,fiE)
      endif
c
      K1=K2
      if(K2.lt.Kmax)goto 100
c
      close(10)
      close(11)
      close(12)
c
      stop
      end
c
!=======================================================================
c
      real*8 FUNCTION fiE(K)
c
      real*8 K,h,x0
      real*8 y,ym1,yp1,x
      real*8 K2
      integer n,Np
c
      external k2
c
      common /param/x0,h
      common /param2/Np
c
      x=x0
      ym1=dsqrt(x0)
      y=dsqrt(x0+h)
c
      do n=1,Np-1
 	 x=x+h
         call Numerow(x,h,y,ym1,yp1,k2,K)
	 ym1=y
	 y=yp1
      enddo
c
      fiE=yp1
c
      return
      end
c
!=======================================================================
c
      SUBROUTINE Numerow(x,h,y,ym1,yp1,k2,K)
c
      implicit none
c
      real*8 x,h,y,ym1,k2,yp1
      real*8 hh,Ap1,Am1,A,K
c
      hh=h**2/12.
      Ap1=1.0+hh*k2(x+h,K)
      Am1=1.0+hh*k2(x-h,K)
      A= 2.0*(1.0-5.0*hh*k2(x,K))
c
      yp1=A/Ap1*y-Am1/Ap1*ym1
c
      return
      end
c
!=======================================================================
c
      real*8 FUNCTION k2(x,K)
c
      implicit none
c
      real*8 x,K
c
      k2=K*K+1./4./x/x
c
      return
      end
c
!=======================================================================
c
      real*8 function ZERO(xl,xr,eps,fun)
c
      implicit none
c
      real*8 xl,xr,xm,eps
      real*8 fun
c
100   xm=(xl+xr)/2.
      if((xr-xl).lt.eps) goto 200
c
      if(fun(xl)*fun(xm).lt.0.0)then
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
